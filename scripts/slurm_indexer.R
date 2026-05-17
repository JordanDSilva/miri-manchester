index_parallel_job_array = function(INDEX, ncores, N, start_Nindex, end_index){
  # will find optimal lots of ncores to fit N galaxies starting at start_Nindex:... inside of a 1 to max_index's slurm JOB array
  # Example N = 200000, ncores = 7, max_index = 800
  # At the minimum I will need floor(20000/(800*7)) = 35 lots of 7 indices 
  # JOB = 1 indices will be 1:7, then 5601:5608, then 11201:11207 ... 35 lots of these < 200000
  # 800*7*35 = 196000 so there will be 4000 remainders 
  # floor(4000/7) = 571 so tack on 7 lots of remainders for JOB = 1-571, and the last 3 stragglers on to JOB=572
  # I could fit 1:100000 galaxies N=100000, then fit 100001:200000 but need to set start_Nindex=100001
  
  maxjobs = end_index * ncores
  how_many_groups = floor(N/maxjobs)
  N_remaining = N - maxjobs*how_many_groups
  
  a = INDEX*ncores - ncores 
  b = INDEX*ncores - 1
  
  if(INDEX == 1){ 
    a = 1 
    b = ncores
    c = 0
  }else{
    c = 1
  }
  job_vec_init = (a:b) + c
  
  job_vec = c()
  for(k in 1:how_many_groups){
    job_vec = c(
      job_vec, job_vec_init + (k-1)*maxjobs
    )
  }
  
  if(N_remaining > 0 ){
    remaining_indices = (maxjobs*how_many_groups + 1):N
    job_vec = as.numeric(na.omit( c(job_vec, remaining_indices[job_vec_init]) )) + (start_Nindex-1)
  }else{
    job_vec = as.numeric(na.omit( c(job_vec) )) + (start_Nindex-1)
  }
  
  return( job_vec ) 
}

index_parallel_job_array(
  1, 
  7,
  14,
  1,
  2
)

N = 200000
ncores = 7
end_index = 800

index_parallel_job_array(
  1, ncores, N = N, end_index = end_index, start_Nindex = 1
)

test = sapply(
  1:end_index,
  function(i){
    index_parallel_job_array(
      i, ncores, N = N, end_index = end_index, start_Nindex = 1
    )
  }
)

indices1 = 1:(N/2)
indices2 = (N/2 + 1):N

ind1 = do.call(
  what = c,
  args = lapply(
    1:end_index,
    function(i){
      index_parallel_job_array(
        i, ncores, N = length(indices1), end_index = end_index, start_Nindex = indices1[1]
      )
    }
  )
)

ind2 = do.call(
  what = c,
  args = lapply(
    1:end_index,
    function(i){
      index_parallel_job_array(
        i, ncores, N = length(indices2), end_index = end_index, start_Nindex = indices2[1]
      )
    }
  )
)

all(
  sort(ind1) == indices1
)
all(
  sort(ind2) == indices2
)
