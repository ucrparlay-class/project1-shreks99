#include <algorithm>
#include "parallel.h"

using namespace parlay;

const size_t THRESHOLD = 10000;
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

//For Calculating prefix sum in parallel
template <class T>
void pscan(T *A, T *B,size_t n) {
  if(n<THRESHOLD) {
    for(size_t i=1;i<n;i++) {
      A[i] = A[i-1] + A[i];
      B[i] = A[i];
    }
    return;
  }
    T* ls = (T*)malloc(n * sizeof(T));
    scan_up(A, ls, n);
    scan_down(A, B, ls, n, 0);
    free(ls);
}

template <class T,class Func>
size_t filter(T* A, size_t n, T* B, size_t pivot,const Func& f) {
    T* flag = (T*)malloc(n * sizeof(T));
    T* ps = (T*)malloc(n * sizeof(T));

    parallel_for(0, n, [&](size_t i) {
        flag[i] = f(A[i],pivot);
    });
    pscan(flag,ps,n);
    parallel_for(0, n, [&](size_t i) {
        if(f(A[i],pivot)) {
          B[ps[i]-1] = A[i]; 
        }
    }
    );
  return ps[n-1];    
}

template <class T>
T para_partition(T *A, size_t n, T pivot) {
  T* B = (T*)malloc(n * sizeof(T)); 
  auto f_less = [](T value, T pivot) { return pivot>value; };
  auto f_more = [](T value, T pivot) { return pivot<value; };
  auto f_equal = [](T value, T pivot) { return pivot==value; };
//less values
  size_t count = filter(A,n,B,pivot,f_less);

//equal values
  count += filter(A,n,B+count,pivot,f_equal);

//large
  filter(A,n,B+count,pivot,f_more);
      
  return count;
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
  if(n<THRESHOLD){
    T pivot = hash64_(std::time(nullptr))%n;
    seq_partition(A,n,pivot);
  }  else
    para_partition(A,n,custom_rand(A,n,THRESHOLD));

}

// template <class T>
// void qsort(T *A, size_t start,size_t end) {
//  // T j;
//   //  if(n<threshold) {
//   //   j = seq_partition(A+start,start,end);
//   // }
//   T  j = para_partition(A+start,end-start+1,A[start+hash64_(std::time(nullptr))%(end-start+1)]);

    
//     par_do(
//         [&] { qsort(A,start, j); },
//        [&] { qsort(A+start,start + j, end); });

// }