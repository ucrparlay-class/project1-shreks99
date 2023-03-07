#include <cmath>
#include <algorithm>
#include "parallel.h"
const size_t THRESHOLD = 10000;

using namespace parlay;

inline uint64_t hash64_(uint64_t u) {
  uint64_t v = u * 3935559000370003845ul + 2691343689449507681ul;
  v ^= v >> 21;
  v ^= v << 37;
  v ^= v >> 4;
  v *= 4768777513237032717ul;
  v ^= v << 20;
  v ^= v >> 41;
  v ^= v << 5;
  return v;
}
template <class T>
void swap(T& a, T& b)
{ 
    T t = a; 
    a = b; 
    b = t; 
}
template <class T>
T seq_partition(T* A, size_t n, T pivot) {
  size_t i = 0, j = n-1;
  while(i<j) {
    while(A[i]<=pivot) i++;
    while(A[j]>pivot) j--;
    if(i<j) swap(A[i],A[j]);
  }
  return j;
}


template <class T>
void scan_down(T* A, T* B, T* ls, size_t n, size_t offset) {
    if(n==1) {
        B[0] = offset + A[0];
        return;
    }
    size_t m = n/2;
    parlay::par_do(
        [&] { scan_down(A,B,ls,m,offset);},
        [&] { scan_down(A + m, B+m, ls + m, n - m,offset + ls[m-1]);});
  }
//For calcuating the leftsum in sequential for prefix-sum
template <class T>
T seq_scan_up(T *A, T *ls, size_t n) {
  T sum = 0;
  for(size_t i=0;i<n;i++)
    sum+=A[i];
  return sum;
}
//For Calulating prefix sum in sequential
template <class T>
void seq_scan(T *A, T *B, size_t n) {
   B[0]=A[0];
    for(size_t i=1;i<n;i++) {
      A[i] = A[i-1] + A[i];
      B[i] = A[i];
    }
    return;
}
//For calculating the leftsum for prefix-sum
template <class T>
T scan_up(T *A, T *ls, size_t n) {
    if(n==1) return A[0];
    size_t m = n/2;
    T l,r;
    auto f1 = [&]() { l = scan_up(A, ls, m);};
    auto f2 = [&]() { r = scan_up(A + m, ls + m, n - m);};
    par_do(f1, f2);

    ls[m-1] = l;
    return l+r;
}
template<typename T>
void pscan_in(T* A, size_t n) {
  if(n==0) return;
  if(n<=THRESHOLD) {
    for(size_t i=1;i<n;i++) A[i] += A[i-1];
    return;
  }
  //From homework 2 prb 5
  // Step 1 Split the array into k chunks, each of size about ⌈n/k⌉ except for probably the last one
  size_t k = ceil(sqrt(n));
  size_t k_size = n/k;
  //cout<<k_size<<" "<<k;
  if(n%k != 0) k_size++;

  //Step 2  sums of k chucks are stored in an array S [1..k]
  T* S = new T[k];
  S[0] = 0;

  //Compute sum of chunks in parallel
  parallel_for(1, k, [&](size_t i) { 
    S[i] = 0;
    for(size_t j = (i-1)*k_size;j<i*k_size;j++) {
      if(j>=n) break;
      S[i] += A[j];
    }
  });
  //Step 3 prefix sum of S[1..k] sequentially
  for(size_t i=1;i<k;i++) S[i] += S[i-1];

  //Step 4 Create k parallel task for k chunks and 
  //compute the prefix sum of a chunk based on input array 
  //and S calculated in previous step
   parallel_for(0, k, [&](size_t i) { 
    for(size_t j = 0;j<k_size;j++) {
      size_t ind = j+k_size*i;
      if(ind < n) {
        if(j==0) A[ind] += S[i];
        else A[ind] += A[ind-1];
      }
    }
  });
  delete[] S;
}

//For Calculating prefix sum in parallel
template <class T>
void pscan(T *A, T *B,size_t n) {
  if(n<THRESHOLD) {
    seq_scan(A,B,n);
  }
    T* ls = (T*)malloc(n * sizeof(T));
    scan_up(A, ls, n);
    scan_down(A, B, ls, n, 0);
    free(ls);
}

template <class T,class Func>
T filter(T* A, size_t n, T* B, T pivot,const Func& f) {
    T* flag = (T*)malloc(n * sizeof(T));
   T* ps = (T*)malloc(n * sizeof(T));
    if(n<THRESHOLD){
      for(size_t i=0;i<n;i++) flag[i] = f(A[i],pivot);
      seq_scan(flag,ps,n);
      for(size_t i=0;i<n;i++) if(f(A[i],pivot)) {
          B[ps[i]-1] = A[i]; 
        }
      return ps[n-1];    
    }
    parallel_for(0, n, [&](size_t i) {
        flag[i] = f(A[i],pivot);
    });
     
   // pscan(flag,ps,n);
    pscan_in(flag,n);
    parallel_for(0, n, [&](size_t i) {
        if(f(A[i],pivot)) {
          B[flag[i]-1] = A[i]; 
        }
    });
    
  return flag[n-1];    
}

template <class T>
T para_partition(T *A, size_t n, T pivot) {
 
  T* B = (T*)malloc(n * sizeof(T)); 
  auto f_less = [](T value, T pivot) { return pivot>value; };
  auto f_more = [](T value, T pivot) { return pivot<value; };
  auto f_equal = [](T value, T pivot) { return pivot==value; };
//less values
  T count = filter(A,n,B,pivot,f_less);
//equal values
  T count2  = count + filter(A,n,B+count,pivot,f_equal);
//large
  filter(A,n,B+count2,pivot,f_more);
  parallel_for(0, n, [&](size_t i) {A[i] = B[i];});
  return count2;
}

template <class T>
T custom_rand(T* A,size_t n,size_t threshold) {
  T* index = (T*)malloc(threshold * sizeof(T));
  for(size_t i=0;i<threshold;i++)
    index[i] = A[hash64_(i)%(n)];
  std::sort(index,index+threshold);
  return index[threshold/2];
}

template <class T>
void quicksort(T *A, size_t n) {
  if(n<=1) return; 
  if(n<THRESHOLD){
    std::sort(A,A+n);
  } else{
  T t = para_partition(A,n,custom_rand(A,n,THRESHOLD));
    if((size_t)t < n){
      parlay::par_do(
            [&] { quicksort(A,t);},
            [&] { quicksort(A+t,n-t);});
    }
  }
}