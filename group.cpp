// Analysis of A_12

#include <other/core/array/Array.h>
#include <other/core/array/convert.h>
#include <other/core/geometry/platonic.h>
#include <other/core/math/factorial.h>
#include <other/core/python/Class.h>
#include <other/core/python/module.h>
#include <other/core/python/numpy.h>
#include <other/core/structure/Hashtable.h>
#include <other/core/utility/Log.h>
#include <other/core/value/Compute.h>
#include <other/core/value/ConstValue.h>
#include <other/core/value/extract.h>
namespace kiss {
namespace {

using namespace other;
using other::to_python;
using Log::cout;
using std::endl;

struct S12 {
  uint64_t p; // little endian mapping from [0,12) to [0,12)

  S12()
    : p(0|1<<4|2<<8|3<<12|4<<16|5<<20|6<<24|7<<28|8l<<32|9l<<36|10l<<40|11l<<44)
  {}

  bool operator==(const S12 x) const {
    return p==x.p;
  }

  int operator[](const int i) const {
    return (p>>4*i)&15;
  }

  S12 operator*(const S12 x) const {
    S12 r;
    r.p = 0;
    for (int i=0;i<12;i++)
      r.p |= uint64_t((*this)[x[i]])<<4*i;
    return r;
  }

  S12 inverse() const {
    S12 r;
    r.p = 0;
    for (int i=0;i<12;i++)
      r.p |= uint64_t(i)<<4*((p>>4*i)&15);
    return r;
  }

  bool valid() const {
    if (p&~((1l<<48)-1))
      return false;
    int mask = 0;
    for (int i=0;i<12;i++)
      mask |= 1l<<((p>>4*i)&15);
    return mask==((1<<12)-1);
  }

  S12 with(const int i, const int pi) const {
    S12 r = *this;
    r.p = (r.p&~(15l<<4*i))|uint64_t(pi)<<4*i;
    return r;
  }
};

PyObject* to_python(const S12 g) {
  return to_python(g.p);
}

S12 from_permutation(RawArray<const int> p) {
  OTHER_ASSERT(p.size()==12);
  S12 r;
  r.p = 0;
  for (const int i : range(12))
    r.p |= uint64_t(p[i])<<4*i;
  OTHER_ASSERT(r.valid());
  return r;
}

Array<const int> to_permutation(const S12 g) {
  Array<int> p(12,false);
  for (const int i : range(12))
    p[i] = g[i];
  return p;
}

int s12_parity(const S12 g) {
  int s = 1;
  int mask = 0;
  for (const int i : range(12))
    if (!(mask&1<<i)) {
      int j = i;
      int n = 0;
      do {
        n++;
        mask |= 1<<j;
        j = g[j]; 
      } while (i != j);
      if (!(n&1))
        s = -s;
    }
  return s;
}

string cycles(const S12 g) {
  string s;
  int mask = 0;
  for (const int i : range(12))
    if (!(mask&1<<i)) {
      if (g[i]!=i) {
        s += '(';
        int j = i;
        do {
          s += j<10?'0'+j:'a'+j-10;
          mask |= 1<<j;
          j = g[j]; 
        } while (i != j);
        s += ')';
      }
    }
  return s;
}

}} namespace other {
using kiss::S12;
template<> struct is_packed_pod<S12> : public mpl::true_ {};
template<> struct NumpyDescr<S12> : public NumpyDescr<uint64_t> {};
template<> struct NumpyIsScalar<S12> : public mpl::true_ {};
ARRAY_CONVERSIONS(1,S12)
template<> struct FromPython<S12> { static S12 convert(PyObject* object) {
  S12 r; 
  r.p = from_python<uint64_t>(object);
  if (!r.valid())
    throw ValueError(format("invalid S12 element %lld",r.p));
  return r;
}};
} namespace kiss { namespace {

S12 s12_times(S12 g, S12 h) { return g*h; }
S12 s12_inverse(S12 g) { return g.inverse(); }

Tuple<Array<const S12>,Array<const Vector<int,2>>> make_generator_info() {
  const auto mesh = icosahedron_mesh().x;
  OTHER_ASSERT(mesh->nodes()==12);
  const auto rings = mesh->sorted_neighbors();
  Array<S12> gens;
  Array<Vector<int,2>> info;
  for (const int i : range(12)) {
    const auto ring = rings[i];
    OTHER_ASSERT(ring.size()==5);
    for (const int n : range(1,5)) {
      S12 g;
      OTHER_ASSERT(g.valid());
      for (const int j : range(5))
        g = g.with(ring[j],ring[(j+n)%5]);
      OTHER_ASSERT(g.valid());
      gens.append(g);
      info.append(vec(i,n<3?n:n-5));
    }
  }
  OTHER_ASSERT(gens.size()==12*4);
  return tuple(gens.const_(),info.const_());
}
const auto generator_info = cache(make_generator_info);
const auto generators = extract<0>(generator_info);

const vector<ValueRef<Hashtable<S12>>> generator_powers = []() {
  vector<ValueRef<Hashtable<S12>>> powers;
  Hashtable<S12> A;
  A.set(S12());
  powers.push_back(const_value(A));
  for (const int k : range(5))
    powers.push_back(cache([=](){
      Hashtable<S12> next;
      const auto gens = generators();
      for (const auto g : powers.at(k)())
        for (const auto h : gens)
          next.set(g*h);    
      return next;
    }));
  return powers;
}();

// Find (k,h) where d(1,g) = k, d(1,h) = k/2, d(h,g) = (k+1)/2
Tuple<int,S12> shortest_path_midpoint(const S12 g) {
  OTHER_ASSERT(s12_parity(g)==1);
  // Check distance 0
  if (g==S12())
    return tuple(0,g);
  // Check distance 1
  if (generator_powers.at(1)().contains(g))
    return tuple(1,S12());
  for (int k=2;;k++) {
    const auto& left = generator_powers.at(k/2)();
    const auto& right = generator_powers.at((k+1)/2)();
    for (const auto h : right)
      if (left.contains(g*h))
        return tuple(k,g*h);
  }
}

Array<S12> shortest_path(const S12 g) {
  const auto kh = shortest_path_midpoint(g);
  Array<S12> path;
  if (kh.x==1)
    path.append(g);
  else if (kh.x>1) {
    path.extend(shortest_path(kh.y));
    path.extend(shortest_path(kh.y.inverse()*g));
  }
  return path;
}

}
}
using namespace kiss;

OTHER_PYTHON_MODULE(kiss_core) {
  OTHER_OBJECT_2(s12_identity,S12())
  OTHER_FUNCTION(s12_parity)
  OTHER_FUNCTION(s12_times)
  OTHER_FUNCTION(s12_inverse)
  OTHER_FUNCTION(from_permutation)
  OTHER_FUNCTION(to_permutation)
  OTHER_FUNCTION(cycles)
  OTHER_OBJECT(generators)
  OTHER_OBJECT(generator_info)
  OTHER_FUNCTION(shortest_path_midpoint)
  OTHER_FUNCTION(shortest_path)
}
