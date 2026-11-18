"""Optional single-device CUDA window aggregation; no GPU RNG or per-base arrays."""
import threading
import time
import numpy as np


def require_cupy():
    try:
        import cupy as cp
        if cp.cuda.runtime.getDeviceCount() < 1:
            raise RuntimeError('No CUDA device')
        return cp
    except (ImportError, RuntimeError) as error:
        raise RuntimeError('GPU backend requires working CuPy/CUDA; choose cpu explicitly if unavailable') from error


_KERNEL = r'''
extern "C" __global__ void windows(
 const long long *ids, const long long *offsets, const long long *ptr,
 const long long *starts, const long long *ends, const long long *first,
 const long long *last, const long long *edges, long long *out, long long *diff,
 const long long count, const long long rows, const long long bins) {
 long long k = (long long)blockIdx.x * blockDim.x + threadIdx.x;
 if(k >= count) return;
 long long lo=0, hi=rows;
 while(lo<hi) { long long m=(lo+hi)/2; if(offsets[m+1]<=k) lo=m+1; else hi=m; }
 long long row=lo, unit=ids[k];
 for(long long b=ptr[unit]; b<ptr[unit+1]; ++b) {
  long long f=first[b], l=last[b];
  if(f==l) atomicAdd((unsigned long long*)&out[row*bins+f], (unsigned long long)(ends[b]-starts[b]));
  else {
   atomicAdd((unsigned long long*)&out[row*bins+f], (unsigned long long)(edges[f+1]-starts[b]));
   atomicAdd((unsigned long long*)&out[row*bins+l], (unsigned long long)(ends[b]-edges[l]));
   atomicAdd((unsigned long long*)&diff[row*(bins+1)+f+1], 1ULL);
   atomicAdd((unsigned long long*)&diff[row*(bins+1)+l], (unsigned long long)-1LL);
  }
 }
}
'''


class GpuWindowAccumulator:
    """Resident block geometry; bounded batch ID transfer and exact int64 atomics.

    CPU sampling threads share one serialized CUDA stream. A batch retains only
    selected unit IDs, batch-by-window counts and difference arrays. No N-by-G
    coverage or dense unit-by-window geometry is built.
    """
    def __init__(self, pool, edges, device=0):
        edges=np.asarray(edges)
        if edges.ndim != 1 or edges.dtype.kind not in 'iu' or len(edges)<2 or edges[0]!=0 or edges[-1]!=pool.genome_len or np.any(edges[1:]<=edges[:-1]):
            raise ValueError('Window edges must partition the full reference')
        self.cp=require_cupy(); self.device=int(device);self.units=len(pool)
        self.bins=len(edges)-1;self.lock=threading.Lock()
        self.calls=0;self.transfer_bytes=0;self.seconds=0.;self.peak_bytes=0
        cp=self.cp
        with cp.cuda.Device(self.device):
            self.starts=cp.asarray(pool.starts,dtype=cp.int64)
            self.ends=cp.asarray(pool.ends,dtype=cp.int64)
            self.ptr=cp.asarray(pool.indptr,dtype=cp.int64)
            self.edges=cp.asarray(edges,dtype=cp.int64)
            self.width=cp.diff(self.edges)
            self.first=cp.searchsorted(self.edges,self.starts,side='right')-1
            self.last=cp.searchsorted(self.edges,self.ends-1,side='right')-1
            self.kernel=cp.RawKernel(_KERNEL,'windows')
            self.kernel.compile()
            cp.cuda.Stream.null.synchronize()
            self.resident_bytes=sum(x.nbytes for x in [self.starts,self.ends,self.ptr,self.edges,self.width,self.first,self.last])
            self.peak_bytes=cp.get_default_memory_pool().used_bytes()

    def count_batch(self, selected):
        return self._count_batch(selected, check_unique=True)

    def _count_selected_batch(self, selected):
        """Internal select_units output already guarantees no repeated IDs.

        Keep cheap shape/range checks and exact mass checks at callers. Public
        count_batch still validates uniqueness for arbitrary caller inputs.
        """
        return self._count_batch(selected, check_unique=False)

    def _count_batch(self, selected, check_unique):
        selected=[np.asarray(ids) for ids in selected]
        if not selected:
            raise ValueError('Require at least one row')
        for ids in selected:
            if ids.ndim!=1 or ids.dtype.kind not in 'iu' or np.any(ids<0) or np.any(ids>=self.units) or (check_unique and len(np.unique(ids))!=len(ids)):
                raise ValueError('Invalid or repeated sampling unit IDs')
        offsets=np.concatenate(([0],np.cumsum([len(x) for x in selected],dtype=np.int64)))
        ids=np.concatenate(selected).astype(np.int64,copy=False)
        cp=self.cp
        with self.lock, cp.cuda.Device(self.device):
            tick=time.perf_counter()
            d_ids=cp.asarray(ids);d_offsets=cp.asarray(offsets)
            out=cp.zeros((len(selected),self.bins),dtype=cp.int64)
            diff=cp.zeros((len(selected),self.bins+1),dtype=cp.int64)
            if len(ids):
                self.kernel(((len(ids)+255)//256,), (256,),
                    (d_ids,d_offsets,self.ptr,self.starts,self.ends,self.first,self.last,self.edges,out,diff,
                     np.int64(len(ids)),np.int64(len(selected)),np.int64(self.bins)))
            out += cp.cumsum(diff[:,:-1],axis=1)*self.width
            self.peak_bytes=max(self.peak_bytes,cp.get_default_memory_pool().used_bytes())
            result=cp.asnumpy(out)  # synchronization and D2H included in elapsed
            self.calls+=1;self.transfer_bytes+=ids.nbytes+offsets.nbytes+result.nbytes
            self.seconds+=time.perf_counter()-tick
            return result

    def metrics(self):
        return dict(device=self.device,cupy=self.cp.__version__,resident_geometry_bytes=self.resident_bytes,
                    peak_live_pool_bytes=self.peak_bytes,reserved_pool_bytes=self.cp.get_default_memory_pool().total_bytes(),
                    transfer_bytes=self.transfer_bytes,calls=self.calls,transfer_and_coverage_seconds=self.seconds,
                    memory_scope='CuPy pool observed after coverage batches, not total device/process high-water mark')
