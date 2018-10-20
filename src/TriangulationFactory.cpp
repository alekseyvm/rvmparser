#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <cassert>
#include "tesselator.h"

#include "Store.h"
#include "Tessellator.h"
#include "GeometryInterface.h"
#include "LinAlgOps.h"

namespace {

  const float pi = float(M_PI);
  const float half_pi = float(0.5*M_PI);
  const float one_over_pi = float(1.0 / M_PI);
  const float twopi = float(2.0*M_PI);

  

 

  // FIXME: replace use of these with stuff from linalg.
  inline void sub3(float* dst, float *a, float *b)
  {
    dst[0] = a[0] - b[0];
    dst[1] = a[1] - b[1];
    dst[2] = a[2] - b[2];
  }

  inline void cross3(float* dst, float *a, float *b)
  {
    dst[0] = a[1] * b[2] - a[2] * b[1];
    dst[1] = a[2] * b[0] - a[0] * b[2];
    dst[2] = a[0] * b[1] - a[1] * b[0];
  }

  inline float dot3(float *a, float *b)
  {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  }

  unsigned triIndices(uint32_t* indices, unsigned  l, unsigned o, unsigned v0, unsigned v1, unsigned v2)
  {
    indices[l++] = o + v0;
    indices[l++] = o + v1;
    indices[l++] = o + v2;
    return l;
  }


  unsigned quadIndices(uint32_t* indices, unsigned  l, unsigned o, unsigned v0, unsigned v1, unsigned v2, unsigned v3)
  {
    indices[l++] = o + v0;
    indices[l++] = o + v1;
    indices[l++] = o + v2;

    indices[l++] = o + v2;
    indices[l++] = o + v3;
    indices[l++] = o + v0;
    return l;
  }

  unsigned quadIndices3(uint32_t* indices, unsigned  t, unsigned o, unsigned v0, unsigned v1, unsigned v2, unsigned v3)
  {
    indices[3 * t + 0] = o + v0;
    indices[3 * t + 1] = o + v1;
    indices[3 * t + 2] = o + v2;

    indices[3 * t + 3] = o + v2;
    indices[3 * t + 4] = o + v3;
    indices[3 * t + 5] = o + v0;
    return t + 2;
  }

  uint32_t* quadIndices_(uint32_t* indices, unsigned o, unsigned v0, unsigned v1, unsigned v2, unsigned v3)
  {
    *indices++ = o + v0;
    *indices++ = o + v1;
    *indices++ = o + v2;

    *indices++ = o + v2;
    *indices++ = o + v3;
    *indices++ = o + v0;

    return indices;
  }


  unsigned vertex(float* normals, float* vertices, unsigned l, float* n, float* p)
  {
    normals[l] = n[0]; vertices[l++] = p[0];
    normals[l] = n[1]; vertices[l++] = p[1];
    normals[l] = n[2]; vertices[l++] = p[2];
    return l;
  }

  unsigned vertex(float* normals, float* vertices, unsigned l, float nx, float ny, float nz, float px, float py, float pz)
  {
    normals[l] = nx; vertices[l++] = px;
    normals[l] = ny; vertices[l++] = py;
    normals[l] = nz; vertices[l++] = pz;
    return l;
  }

  unsigned tessellateCircle(uint32_t* indices, unsigned  l, uint32_t* t, uint32_t* src, unsigned N)
  {
    while (3 <= N) {
      unsigned m = 0;
      unsigned i;
      for (i = 0; i + 2 < N; i += 2) {
        indices[l++] = src[i];
        indices[l++] = src[i + 1];
        indices[l++] = src[i + 2];
        t[m++] = src[i];
      }
      for (; i < N; i++) {
        t[m++] = src[i];
      }
      N = m;
      std::swap(t, src);
    }
    return l;
  }

  unsigned tessellateCircle3(uint32_t* indices, unsigned  t, uint32_t* tmp, uint32_t* src, unsigned N)
  {
    while (3 <= N) {
      unsigned m = 0;
      unsigned i;
      for (i = 0; i + 2 < N; i += 2) {
        indices[3 * t + 0] = src[i];
        indices[3 * t + 1] = src[i + 1];
        indices[3 * t + 2] = src[i + 2];
        t++;
        tmp[m++] = src[i];
      }
      for (; i < N; i++) {
        tmp[m++] = src[i];
      }
      N = m;
      std::swap(tmp, src);
    }
    return t;
  }

  
  uint32_t* tessellateCircle_(uint32_t* I, unsigned offset, unsigned flip,  unsigned N)
  {
    auto I0 = I;
    for (unsigned step = 2; step < N; step = 2 * step) {
      for (unsigned i = 0; i + step / 2 < N; i += step) {
        *I++ = offset + i + (flip ? 0 : step / 2);
        *I++ = offset + i + (flip ? step / 2 : 0);
        *I++ = offset + (i + step < N ? i + step : 0);
      }
    }
    return I;
  }

  Vec2f encodeCylindricalTexCoord(const float unitCircum, const float distance, const float circumference)
  {
    return Vec2f(unitCircum, distance);
  }

  // tex in circle on origin with radius 1.
  Vec2f encodeCircularTexCoord(const Vec2f& tex, const float radius)
  {
    return radius * tex;
  }


}


TriangulationFactory::TriangulationFactory(Store* store, Logger logger, float tolerance, unsigned minSamples, unsigned maxSamples, bool texcoords, bool smoothingGroups) :
  store(store),
  logger(logger),
  tolerance(tolerance),
  minSamples(minSamples),
  maxSamples(std::max(minSamples, maxSamples)),
  texcoords(texcoords),
  smoothingGroups(smoothingGroups),
  nextSmoothingGroup(store->componentCount())
{
}


unsigned TriangulationFactory::sagittaBasedSegmentCount(float arc, float radius, float scale)
{
  float samples = arc / std::acos(std::max(-1.f, 1.f - tolerance / (scale*radius)));
  return std::min(maxSamples, unsigned(std::max(float(minSamples), std::ceil(samples))));
}


float TriangulationFactory::sagittaBasedError(float arc, float radius, float scale, unsigned segments)
{
  auto s = scale * radius*(1.f - std::cos(arc / segments));  // Length of sagitta
  //assert(s <= tolerance);
  return s;
}




Triangulation* TriangulationFactory::pyramid(Arena* arena, const Geometry* geo, float scale)
{
  auto bx = 0.5f * geo->pyramid.bottom[0];
  auto by = 0.5f * geo->pyramid.bottom[1];
  auto tx = 0.5f * geo->pyramid.top[0];
  auto ty = 0.5f * geo->pyramid.top[1];
  auto ox = 0.5f * geo->pyramid.offset[0];
  auto oy = 0.5f * geo->pyramid.offset[1];
  auto h2 = 0.5f * geo->pyramid.height;

  

  Vec3f quad[2][4] =
  {
    {
      Vec3f( -bx - ox, -by - oy, -h2 ),
      Vec3f(  bx - ox, -by - oy, -h2 ),
      Vec3f(  bx - ox,  by - oy, -h2 ),
      Vec3f(-bx - ox,  by - oy, -h2 )
    },
    {
       Vec3f(-tx + ox, -ty + oy, h2),
       Vec3f(tx + ox, -ty + oy, h2),
       Vec3f(tx + ox,  ty + oy, h2),
       Vec3f(-tx + ox,  ty + oy, h2)
    },
  };

  Vec3f n[6] = {
    Vec3f( 0.f, -h2,  (quad[1][0][1] - quad[0][0][1]) ),
    Vec3f(  h2, 0.f, -(quad[1][1][0] - quad[0][1][0]) ),
    Vec3f( 0.f,  h2, -(quad[1][2][1] - quad[0][2][1]) ),
    Vec3f( -h2, 0.f,  (quad[1][3][0] - quad[0][3][0]) ),
    Vec3f(0, 0, -1 ),
    Vec3f(0, 0, 1),
  };

  bool cap[6] = {
      true,
      true,
      true,
      true,
      1e-7f <= std::min(std::abs(geo->pyramid.bottom[0]), std::abs(geo->pyramid.bottom[1])),
      1e-7f <= std::min(std::abs(geo->pyramid.top[0]), std::abs(geo->pyramid.top[1]))
  };

  for (unsigned i = 0; i < 6; i++) {
    auto * con = geo->connections[i];
    if (cap[i] == false || con == nullptr || con->flags != Connection::Flags::HasRectangularSide) continue;

    if (doInterfacesMatch(geo, con)) {
      cap[i] = false;
      discardedCaps++;
      //store->addDebugLine(con->p.data, (con->p.data + 0.05f*con->d).data, 0xff0000);
    }
  }

  unsigned caps = 0;
  for (unsigned i = 0; i < 6; i++) if (cap[i]) caps++;

  auto * tri = arena->alloc<Triangulation>();
  tri->error = 0.f;

  tri->vertices_n = 4 * caps;
  tri->triangles_n = 2 * caps;

  tri->vertices = (float*)arena->alloc(3 * sizeof(float)*tri->vertices_n);
  tri->normals = (float*)arena->alloc(3 * sizeof(float)*tri->vertices_n);

  unsigned l = 0;
  for (unsigned i = 0; i < 4; i++) {
    if (cap[i] == false) continue;
    unsigned ii = (i + 1) & 3;
    l = vertex(tri->normals, tri->vertices, l, n[i].data, quad[0][i].data);
    l = vertex(tri->normals, tri->vertices, l, n[i].data, quad[0][ii].data);
    l = vertex(tri->normals, tri->vertices, l, n[i].data, quad[1][ii].data);
    l = vertex(tri->normals, tri->vertices, l, n[i].data, quad[1][i].data);
  }
  if (cap[4]) {
    for (unsigned i = 0; i < 4; i++) {
      l = vertex(tri->normals, tri->vertices, l, n[4].data, quad[0][i].data);
    }
  }
  if (cap[5]) {
    for (unsigned i = 0; i < 4; i++) {
      l = vertex(tri->normals, tri->vertices, l, n[5].data, quad[1][i].data);
    }
  }
  assert(l == 3*tri->vertices_n);

  l = 0;
  unsigned o = 0;
  tri->indices = (uint32_t*)arena->alloc(3 * sizeof(uint32_t) * tri->triangles_n);
  for (unsigned i = 0; i < 4; i++) {
    if (cap[i] == false) continue;
    l = quadIndices(tri->indices, l, o /*4 * i*/, 0, 1, 2, 3);
    o += 4;
  }
  if (cap[4]) {
    l = quadIndices(tri->indices, l, o, 3, 2, 1, 0);
    o += 4;
  }
  if (cap[5]) {
    l = quadIndices(tri->indices, l, o, 0, 1, 2, 3);
    o += 4;
  }
  assert(l == 3 * tri->triangles_n);
  assert(o == tri->vertices_n);

  return tri;
}


Triangulation* TriangulationFactory::box(Arena* arena, const Geometry* geo, float scale)
{
  auto & box = geo->box;

  auto xp = 0.5f * box.lengths[0]; auto xm = -xp;
  auto yp = 0.5f * box.lengths[1]; auto ym = -yp;
  auto zp = 0.5f * box.lengths[2]; auto zm = -zp;

  Vec3f V[6][4] = {
    { Vec3f(xm, ym, zp ), Vec3f(xm, yp, zp ), Vec3f(xm, yp, zm ), Vec3f(xm, ym, zm ) },
    { Vec3f(xp, ym, zm ), Vec3f(xp, yp, zm ), Vec3f(xp, yp, zp ), Vec3f(xp, ym, zp ) },
    { Vec3f(xp, ym, zm ), Vec3f(xp, ym, zp ), Vec3f(xm, ym, zp ), Vec3f(xm, ym, zm ) },
    { Vec3f(xm, yp, zm ), Vec3f(xm, yp, zp ), Vec3f(xp, yp, zp ), Vec3f(xp, yp, zm ) },
    { Vec3f(xm, yp, zm ), Vec3f(xp, yp, zm ), Vec3f(xp, ym, zm ), Vec3f(xm, ym, zm ) },
    { Vec3f(xm, ym, zp ), Vec3f(xp, ym, zp ), Vec3f(xp, yp, zp ), Vec3f(xm, yp, zp ) }
  };

  Vec3f N[6] = {
    Vec3f(-1,  0,  0),
    Vec3f(1,  0,  0 ),
    Vec3f(0, -1,  0 ),
    Vec3f(0,  1,  0 ),
    Vec3f(0,  0, -1 ),
    Vec3f(0,  0,  1 )
  };

  bool faces[6] = {
    1e-5 <= box.lengths[0],
    1e-5 <= box.lengths[0],
    1e-5 <= box.lengths[1],
    1e-5 <= box.lengths[1],
    1e-5 <= box.lengths[2],
    1e-5 <= box.lengths[2],
  };
  for (unsigned i = 0; i < 6; i++) {
    auto * con = geo->connections[i];
    if (faces[i] == false || con == nullptr || con->flags != Connection::Flags::HasRectangularSide) continue;

    if (doInterfacesMatch(geo, con)) {
      faces[i] = false;
      discardedCaps++;
      //store->addDebugLine(con->p.data, (con->p.data + 0.05f*con->d).data, 0xff0000);
    }

  }

  unsigned faces_n = 0;
  for(unsigned i=0; i<6; i++) {
    if (faces[i]) faces_n++;
  }

  Triangulation* tri = arena->alloc<Triangulation>();
  tri->error = 0.f;

  if (faces_n) {
    tri->vertices_n = 4 * faces_n;
    tri->vertices = (float*)arena->alloc(3 * sizeof(float)*tri->vertices_n);
    tri->normals = (float*)arena->alloc(3 * sizeof(float)*tri->vertices_n);

    tri->triangles_n = 2 * faces_n;
    tri->indices = (uint32_t*)arena->alloc(3 * sizeof(uint32_t)*tri->triangles_n);

    unsigned o = 0;
    unsigned i_v = 0;
    unsigned i_p = 0;
    for (unsigned f = 0; f < 6; f++) {
      if (!faces[f]) continue;

      for (unsigned i = 0; i < 4; i++) {
        i_v = vertex(tri->normals, tri->vertices, i_v, N[f].data, V[f][i].data);
      }
      i_p = quadIndices(tri->indices, i_p, o, 0, 1, 2, 3);

      o += 4;
    }
    assert(i_v == 3 * tri->vertices_n);
    assert(i_p == 3 * tri->triangles_n);
    assert(o == tri->vertices_n);
  }

  return tri;
}


Triangulation* TriangulationFactory::rectangularTorus(Arena* arena, const Geometry* geo, float scale)
{
  //if (cullTiny && std::max(tor.outer_radius-tor.inner_radius, tor.height)*scale < tolerance) {
  //  tri->error = std::max(tor.outer_radius - tor.inner_radius, tor.height)*scale;
  //  return;
  //}

  auto & tor = geo->rectangularTorus;
  auto segments = sagittaBasedSegmentCount(tor.angle, tor.outer_radius, scale);
  auto samples = segments + 1;  // Assumed to be open, add extra sample.

  auto * tri = arena->alloc<Triangulation>();
  tri->error = sagittaBasedError(tor.angle, tor.outer_radius, scale, segments);

  bool shell = true;
  bool cap[2] = {
    true,
    true
  };

  for (unsigned i = 0; i < 2; i++) {
    auto * con = geo->connections[i];
    if (con && con->flags == Connection::Flags::HasRectangularSide) {
      if (doInterfacesMatch(geo, con)) {
        cap[i] = false;
        discardedCaps++;
        //store->addDebugLine(con->p.data, (con->p.data + 0.05f*con->d).data, 0xff0000);
      }
    }
  }

  auto h2 = 0.5f*tor.height;
  float square[4][2] = {
    { tor.outer_radius, -h2 },
    { tor.inner_radius, -h2 },
    { tor.inner_radius, h2 },
    { tor.outer_radius, h2 },
  };

  // Not closed
  t0.resize(2 * samples + 1);
  for (unsigned i = 0; i < samples; i++) {
    t0[2 * i + 0] = std::cos((tor.angle / segments)*i);
    t0[2 * i + 1] = std::sin((tor.angle / segments)*i);
  }

  unsigned l = 0;

  tri->vertices_n = (shell ? 4 * 2 * samples : 0) + (cap[0] ? 4 : 0) + (cap[1] ? 4 : 0);
  tri->vertices = (float*)arena->alloc(3 * sizeof(float)*tri->vertices_n);
  tri->normals = (float*)arena->alloc(3 * sizeof(float)*tri->vertices_n);

  if (shell) {
    for (unsigned i = 0; i < samples; i++) {
      float n[4][3] = {
        { 0.f, 0.f, -1.f },
        { -t0[2 * i + 0], -t0[2 * i + 1], 0.f },
        { 0.f, 0.f, 1.f },
        { t0[2 * i + 0], t0[2 * i + 1], 0.f },
      };

      for (unsigned k = 0; k < 4; k++) {
        unsigned kk = (k + 1) & 3;

        tri->normals[l] = n[k][0]; tri->vertices[l++] = square[k][0] * t0[2 * i + 0];
        tri->normals[l] = n[k][1]; tri->vertices[l++] = square[k][0] * t0[2 * i + 1];
        tri->normals[l] = n[k][2]; tri->vertices[l++] = square[k][1];

        tri->normals[l] = n[k][0]; tri->vertices[l++] = square[kk][0] * t0[2 * i + 0];
        tri->normals[l] = n[k][1]; tri->vertices[l++] = square[kk][0] * t0[2 * i + 1];
        tri->normals[l] = n[k][2]; tri->vertices[l++] = square[kk][1];
      }
    }
  }
  if (cap[0]) {
    for (unsigned k = 0; k < 4; k++) {
      tri->normals[l] =  0.f; tri->vertices[l++] = square[k][0] * t0[0];
      tri->normals[l] = -1.f; tri->vertices[l++] = square[k][0] * t0[1];
      tri->normals[l] =  0.f; tri->vertices[l++] = square[k][1];
    }
  }
  if (cap[1]) {
    for (unsigned k = 0; k < 4; k++) {
      tri->normals[l] = -t0[2 * (samples - 1) + 1]; tri->vertices[l++] = square[k][0] * t0[2 * (samples - 1) + 0];
      tri->normals[l] =  t0[2 * (samples - 1) + 0]; tri->vertices[l++] = square[k][0] * t0[2 * (samples - 1) + 1];
      tri->normals[l] =                        0.f; tri->vertices[l++] = square[k][1];
    }
  }
  assert(l == 3*tri->vertices_n);

  l = 0;
  unsigned o = 0;

  tri->triangles_n = (shell ? 4 * 2 * (samples - 1) : 0) + (cap[0] ? 2 : 0) + (cap[1] ? 2 : 0);
  tri->indices = (uint32_t*)arena->alloc(3 * sizeof(uint32_t)*tri->triangles_n);

  if (shell) {
    for (unsigned i = 0; i + 1 < samples; i++) {
      for (unsigned k = 0; k < 4; k++) {
        tri->indices[l++] = 4 * 2 * (i + 0) + 0 + 2 * k;
        tri->indices[l++] = 4 * 2 * (i + 0) + 1 + 2 * k;
        tri->indices[l++] = 4 * 2 * (i + 1) + 0 + 2 * k;

        tri->indices[l++] = 4 * 2 * (i + 1) + 0 + 2 * k;
        tri->indices[l++] = 4 * 2 * (i + 0) + 1 + 2 * k;
        tri->indices[l++] = 4 * 2 * (i + 1) + 1 + 2 * k;
      }
    }
    o += 4 * 2 * samples;
  }
  if (cap[0]) {
    tri->indices[l++] = o + 0;
    tri->indices[l++] = o + 2;
    tri->indices[l++] = o + 1;
    tri->indices[l++] = o + 2;
    tri->indices[l++] = o + 0;
    tri->indices[l++] = o + 3;
    o += 4;
  }
  if (cap[1]) {
    tri->indices[l++] = o + 0;
    tri->indices[l++] = o + 1;
    tri->indices[l++] = o + 2;
    tri->indices[l++] = o + 2;
    tri->indices[l++] = o + 3;
    tri->indices[l++] = o + 0;
    o += 4;
  }
  assert(o == tri->vertices_n);
  assert(l == 3*tri->triangles_n);

  return tri;
}


Triangulation* TriangulationFactory::circularTorus(Arena* arena, const Geometry* geo, float scale)
{
  auto & ct = geo->circularTorus;
  unsigned segments_l = sagittaBasedSegmentCount(ct.angle, ct.offset + ct.radius, scale); // large radius, toroidal direction
  unsigned segments_s = sagittaBasedSegmentCount(twopi, ct.radius, scale); // small radius, poloidal direction

  auto * tri = arena->alloc<Triangulation>();
  tri->error = std::max(sagittaBasedError(ct.angle, ct.offset + ct.radius, scale, segments_l),
                        sagittaBasedError(twopi, ct.radius, scale, segments_s));

  unsigned samplesMajor = segments_l + 1;  // Assumed to be open, add extra sample
  unsigned samplesMinor = segments_s;      // Assumed to be closed

  unsigned shell = 1;
  bool cap[2] = { true, true };
  unsigned caps = 0;
  for (unsigned i = 0; i < 2; i++) {
    auto * con = geo->connections[i];
    if (con && con->flags == Connection::Flags::HasCircularSide) {
      if (doInterfacesMatch(geo, con)) {
        cap[i] = false;
        discardedCaps++;
      }
    }
    if (cap[i]) caps++;
  }

  cosSinMajor.resize(samplesMajor); // t0
  for (unsigned i = 0; i < samplesMajor; i++) {
    cosSinMajor[i] = Vec2f(std::cos((ct.angle / (samplesMajor - 1.f))*i),
                           std::sin((ct.angle / (samplesMajor - 1.f))*i));
  }

  cosSin.resize(samplesMinor);  // t1
  for (unsigned i = 0; i < samplesMinor; i++) {
    cosSin[i] = Vec2f(std::cos((twopi / samplesMinor)*i + geo->sampleStartAngle),
                      std::sin((twopi / samplesMinor)*i + geo->sampleStartAngle));
  }


  if (texcoords == false) {
    tri->vertices_n = samplesMinor * (samplesMajor * shell + caps);
    tri->texCoords = nullptr;
  }
  else {
    tri->vertices_n = (samplesMinor + 1)*samplesMajor * shell + samplesMinor * caps;
    tri->texCoords = (float*)arena->alloc(sizeof(Vec2f)*tri->vertices_n);
  }

  tri->vertices = (float*)arena->alloc(3 * sizeof(float)*tri->vertices_n);
  tri->normals = (float*)arena->alloc(3 * sizeof(float)*tri->vertices_n);
  tri->triangles_n = 2 * (samplesMajor - 1)*samplesMinor * shell + (samplesMinor - 2)*caps;
  tri->indices = (uint32_t*)arena->alloc(3 * sizeof(uint32_t)*tri->triangles_n);
  if (smoothingGroups) {
    tri->smoothingGroups = (uint32_t*)arena->alloc(sizeof(uint32_t)*tri->triangles_n);
  }

  auto * V = (Vec3f*)tri->vertices;
  auto * N = (Vec3f*)tri->normals;
  auto * T = (Vec2f*)tri->texCoords;
  auto * I = tri->indices;
  auto * S = tri->smoothingGroups;
  unsigned o = 0;

  if (texcoords == false) {
    if (shell) {
      for (unsigned u = 0; u < samplesMajor; u++) {
        for (unsigned v = 0; v < samplesMinor; v++) {
          *V++ = Vec3f((ct.radius * cosSin[v].x + ct.offset) * cosSinMajor[u], ct.radius * cosSin[v].y);
          *N++ = Vec3f(cosSin[v].x * cosSinMajor[u], cosSin[v].y);
        }
      }
    }
    if (cap[0]) {
      for (unsigned v = 0; v < samplesMinor; v++) {
        *V++ = Vec3f((ct.radius*cosSin[v].x + ct.offset) * cosSinMajor[0], ct.radius*cosSin[v].y);
        *N++ = Vec3f(0, -1, 0);
      }
    }
    if (cap[1]) {
      unsigned m = samplesMajor - 1;
      for (unsigned v = 0; v < samplesMinor; v++) {
        *V++ = Vec3f((ct.radius * cosSin[v].x + ct.offset) * cosSinMajor[m], ct.radius * cosSin[v].y);
        *N++ = Vec3f(-cosSinMajor[m].y, cosSinMajor[m].x, 0.f);
      }
    }
    if (shell) {
      for (unsigned u = 0; u + 1 < samplesMajor; u++) {
        for (unsigned v = 0; v + 1 < samplesMinor; v++) {
          I = quadIndices_(I, o,
                           samplesMinor * (u + 0) + (v + 0),
                           samplesMinor * (u + 1) + (v + 0),
                           samplesMinor * (u + 1) + (v + 1),
                           samplesMinor * (u + 0) + (v + 1));
        }
        I = quadIndices_(I, o,
                         samplesMinor * (u + 0) + (samplesMinor - 1),
                         samplesMinor * (u + 1) + (samplesMinor - 1),
                         samplesMinor * (u + 1) + 0,
                         samplesMinor * (u + 0) + 0);
      }
      o += samplesMajor * samplesMinor;
      if (smoothingGroups) {
        for (unsigned i = 0; i < 2 * (samplesMajor - 1)*samplesMinor; i++) {
          *S++ = geo->componentId;
        }
      }
    }
    if (cap[0]) {
      I = tessellateCircle_(I, o, 1, samplesMinor);
      o += samplesMinor;
      if (smoothingGroups) {
        auto sgrp = newSmoothingGroup();
        for (unsigned i = 0; i < (samplesMinor - 2); i++) {
          *S++ = sgrp;
        }
      }
    }
    if (cap[1]) {
      I = tessellateCircle_(I, o, 0, samplesMinor);
      o += samplesMinor;
      if (smoothingGroups) {
        auto sgrp = newSmoothingGroup();
        for (unsigned i = 0; i < (samplesMinor - 2); i++) {
          *S++ = sgrp;
        }
      }
    }
  }
  else {
    const float distanceA = geo->distances[0];
    const float distanceB = geo->distances[1];
    const auto circumference = scale * twopi * ct.radius;
    if (shell) {
      float distances[2] = {
        geo->distances[0],
        geo->distances[1]
      };
      const auto vScale = distances[0] < distances[1] ? -1.f / samplesMinor : 1.f / samplesMinor;
      const auto vConst = distances[0] < distances[1] ? 1 : 0;
      const auto uScale = 1.f / (samplesMajor - 1);
      for (unsigned u = 0; u < samplesMajor; u++) {
        const auto un = uScale * u;
        const auto distance = (1.f - un)*distanceA + un * distanceB;
        for (unsigned v = 0; v < samplesMinor; v++) {
          *V++ = Vec3f((ct.radius * cosSin[v].x + ct.offset) * cosSinMajor[u], ct.radius * cosSin[v].y);
          *N++ = Vec3f(cosSin[v].x * cosSinMajor[u], cosSin[v].y);
          *T++ = encodeCylindricalTexCoord(vScale * v + vConst, distance, circumference);
        }
        *V++ = Vec3f((ct.radius * cosSin[0].x + ct.offset) * cosSinMajor[u], ct.radius * cosSin[0].y);
        *N++ = Vec3f(cosSin[0].x * cosSinMajor[u], cosSin[0].y);
        *T++ = encodeCylindricalTexCoord(1.f - vConst, distance, circumference);
      }
    }
    if (cap[0]) {
      for (unsigned v = 0; v < samplesMinor; v++) {
        *V++ = Vec3f((ct.radius*cosSin[v].x + ct.offset) * cosSinMajor[0], ct.radius*cosSin[v].y);
        *N++ = Vec3f(0, -1, 0);
        *T++ = encodeCircularTexCoord(cosSin[v], scale * ct.radius);
      }
    }
    if (cap[1]) {
      unsigned m = samplesMajor - 1;
      for (unsigned v = 0; v < samplesMinor; v++) {
        *V++ = Vec3f((ct.radius * cosSin[v].x + ct.offset) * cosSinMajor[m], ct.radius * cosSin[v].y);
        *N++ = Vec3f(-cosSinMajor[m].y, cosSinMajor[m].x, 0.f);
        *T++ = encodeCircularTexCoord(cosSin[v], scale * ct.radius);
      }
    }
    if (shell) {
      for (unsigned u = 0; u + 1 < samplesMajor; u++) {
        for (unsigned v = 0; v < samplesMinor; v++) {
          I = quadIndices_(I, o,
                           (samplesMinor + 1) * (u + 0) + (v + 0),
                           (samplesMinor + 1) * (u + 1) + (v + 0),
                           (samplesMinor + 1) * (u + 1) + (v + 1),
                           (samplesMinor + 1) * (u + 0) + (v + 1));
        }
      }
      o += samplesMajor * (samplesMinor + 1);
      if (smoothingGroups) {
        for (unsigned i = 0; i < 2 * (samplesMajor - 1)*samplesMinor; i++) {
          *S++ = geo->componentId;
        }
      }
    }
    if (cap[0]) {
      I = tessellateCircle_(I, o, 1, samplesMinor);
      o += samplesMinor;
      if (smoothingGroups) {
        auto sgrp = newSmoothingGroup();
        for (unsigned i = 0; i < (samplesMinor - 2); i++) {
          *S++ = sgrp;
        }
      }
    }
    if (cap[1]) {
      I = tessellateCircle_(I, o, 0, samplesMinor);
      o += samplesMinor;
      if (smoothingGroups) {
        auto sgrp = newSmoothingGroup();
        for (unsigned i = 0; i < (samplesMinor - 2); i++) {
          *S++ = sgrp;
        }
      }
    }

  }

  assert((V - (Vec3f*)tri->vertices) == tri->vertices_n);
  assert((N - (Vec3f*)tri->normals) == tri->vertices_n);
  assert((I - tri->indices) == 3 * tri->triangles_n);
  assert(o == tri->vertices_n);
  if (smoothingGroups) {
    assert((S - tri->smoothingGroups) == tri->triangles_n);
  }
  return tri;
}


Triangulation* TriangulationFactory::snout(Arena* arena, const Geometry* geo, float scale)
{
  //if (cullTiny && radius_max*scale < tolerance) {
  //  tri->error = radius_max *scale;
  //  return;
  //}
  
  auto & sn = geo->snout;
  auto radius_max = std::max(sn.radius_b, sn.radius_t);
  unsigned segments = sagittaBasedSegmentCount(twopi, radius_max, scale);
  unsigned samples = segments;  // assumed to be closed

  auto * tri = arena->alloc<Triangulation>();
  tri->error = sagittaBasedError(twopi, radius_max, scale, segments);

  bool shell = true;
  bool cap[2] = { true, true };
  float radii[2] = { geo->snout.radius_b, geo->snout.radius_t };
  for (unsigned i = 0; i < 2; i++) {
    auto * con = geo->connections[i];
    if (con && con->flags == Connection::Flags::HasCircularSide) {
      if (doInterfacesMatch(geo, con)) {
        cap[i] = false;
        discardedCaps++;
      }
      else {
        //store->addDebugLine(con->p.data, (con->p.data + 0.05f*con->d).data, 0x00ffff);
      }
    }
  }

  t0.resize(2 * samples);
  for (unsigned i = 0; i < samples; i++) {
    t0[2 * i + 0] = std::cos((twopi / samples)*i + geo->sampleStartAngle);
    t0[2 * i + 1] = std::sin((twopi / samples)*i + geo->sampleStartAngle);
  }
  t1.resize(2 * samples);
  for (unsigned i = 0; i < 2 * samples; i++) {
    t1[i] = sn.radius_b * t0[i];
  }
  t2.resize(2 * samples);
  for (unsigned i = 0; i < 2 * samples; i++) {
    t2[i] = sn.radius_t * t0[i];
  }

  float h2 = 0.5f*sn.height;
  unsigned l = 0;
  auto ox = 0.5f*sn.offset[0];
  auto oy = 0.5f*sn.offset[1];
  float mb[2] = { std::tan(sn.bshear[0]), std::tan(sn.bshear[1]) };
  float mt[2] = { std::tan(sn.tshear[0]), std::tan(sn.tshear[1]) };

  tri->vertices_n = (shell ? 2 * samples : 0) + (cap[0] ? samples : 0) + (cap[1] ? samples : 0);
  tri->vertices = (float*)arena->alloc(3 * sizeof(float)*tri->vertices_n);
  tri->normals = (float*)arena->alloc(3 * sizeof(float)*tri->vertices_n);

  tri->triangles_n = (shell ? 2 * samples : 0) + (cap[0] ? samples - 2 : 0) + (cap[1] ? samples - 2 : 0);
  tri->indices = (uint32_t*)arena->alloc(3 * sizeof(uint32_t)*tri->triangles_n);

  if (shell) {
    for (unsigned i = 0; i < samples; i++) {
      float xb = t1[2 * i + 0] - ox;
      float yb = t1[2 * i + 1] - oy;
      float zb = -h2 + mb[0] * t1[2 * i + 0] + mb[1] * t1[2 * i + 1];

      float xt = t2[2 * i + 0] + ox;
      float yt = t2[2 * i + 1] + oy;
      float zt = h2 + mt[0] * t2[2 * i + 0] + mt[1] * t2[2 * i + 1];

      float s = (sn.offset[0] * t0[2 * i + 0] + sn.offset[1] * t0[2 * i + 1]);
      float nx = t0[2 * i + 0];
      float ny = t0[2 * i + 1];
      float nz = -(sn.radius_t - sn.radius_b + s) / sn.height;

      l = vertex(tri->normals, tri->vertices, l, nx, ny, nz, xb, yb, zb);
      l = vertex(tri->normals, tri->vertices, l, nx, ny, nz, xt, yt, zt);
    }
  }
  if (cap[0]) {
    auto nx = std::sin(sn.bshear[0])*std::cos(sn.bshear[1]);
    auto ny = std::sin(sn.bshear[1]);
    auto nz = -std::cos(sn.bshear[0])*std::cos(sn.bshear[1]);
    for (unsigned i = 0; cap[0] && i < samples; i++) {
      l = vertex(tri->normals, tri->vertices, l, nx, ny, nz,
                 t1[2 * i + 0] - ox,
                 t1[2 * i + 1] - oy,
                 -h2 + mb[0] * t1[2 * i + 0] + mb[1] * t1[2 * i + 1]);
    }
  }
  if (cap[1]) {
    auto nx = -std::sin(sn.tshear[0])*std::cos(sn.tshear[1]);
    auto ny = -std::sin(sn.tshear[1]);
    auto nz = std::cos(sn.tshear[0])*std::cos(sn.tshear[1]);
    for (unsigned i = 0; i < samples; i++) {
      l = vertex(tri->normals, tri->vertices, l, nx, ny, nz,
                 t2[2 * i + 0] + ox,
                 t2[2 * i + 1] + oy,
                 h2 + mt[0] * t2[2 * i + 0] + mt[1] * t2[2 * i + 1]);
    }
  }
  assert(l == 3 * tri->vertices_n);

  l = 0;
  unsigned o = 0;
  if (shell) {
    for (unsigned i = 0; i < samples; i++) {
      unsigned ii = (i + 1) % samples;
      l = quadIndices(tri->indices, l, 0, 2 * i, 2 * ii, 2 * ii + 1, 2 * i + 1);
    }
    o += 2 * samples;
  }

  u1.resize(samples);
  u2.resize(samples);
  if (cap[0]) {
    for (unsigned i = 0; i < samples; i++) {
      u1[i] = o + (samples - 1) - i;
    }
    l = tessellateCircle(tri->indices, l, u2.data(), u1.data(), samples);
    o += samples;
  }
  if (cap[1]) {
    for (unsigned i = 0; i < samples; i++) {
      u1[i] = o + i;
    }
    l = tessellateCircle(tri->indices, l, u2.data(), u1.data(), samples);
    o += samples;
  }
  assert(l == 3*tri->triangles_n);
  assert(o == tri->vertices_n);

  return tri;
}


Triangulation* TriangulationFactory::cylinder(Arena* arena, const Geometry* geo, float scale)
{
  const auto & cy = geo->cylinder;
  unsigned segments = sagittaBasedSegmentCount(twopi, cy.radius, scale);
  unsigned samples = segments;  // Assumed to be closed


  auto * tri = arena->alloc<Triangulation>();
  tri->error = sagittaBasedError(twopi, cy.radius, scale, segments);

  unsigned shell = 1;
  bool cap[2] = { true, true };
  unsigned caps = 0;
  for (unsigned i = 0; i < 2; i++) {
    auto * con = geo->connections[i];
    if (con && con->flags == Connection::Flags::HasCircularSide) {
      if (doInterfacesMatch(geo, con)) {
        cap[i] = false;
        discardedCaps++;
      }
    }
    if (cap[i]) caps++;
  }

  cosSin.resize(samples);
  for (unsigned i = 0; i < samples; i++) {
    cosSin[i] = Vec2f(std::cos((twopi / samples)*i + geo->sampleStartAngle),
                      std::sin((twopi / samples)*i + geo->sampleStartAngle));
  }
  cosSinRadius.resize(samples);
  for (unsigned i = 0; i < samples; i++) {
    cosSinRadius[i] = cy.radius * cosSin[i];
  }

  if (texcoords == false) {
    tri->vertices_n = samples * ((shell ? 2 : 0) + (cap[0] ? 1 : 0) + (cap[1] ? 1 : 0));
    tri->triangles_n = (shell ? 2 * samples : 0) + (cap[0] ? samples - 2 : 0) + (cap[1] ? samples - 2 : 0);
    tri->texCoords = nullptr;
  }
  else {
    tri->vertices_n = 2 * (samples + 1) * shell + samples * caps;
    tri->triangles_n = 2 * samples*shell + (samples - 2)*caps;
    tri->texCoords = (float*)arena->alloc(sizeof(Vec2f)*tri->vertices_n);
  }
  if (smoothingGroups) {
    tri->smoothingGroups = (uint32_t*)arena->alloc(sizeof(uint32_t)*tri->triangles_n);
  }

  tri->vertices = (float*)arena->alloc(sizeof(Vec3f)*tri->vertices_n);
  tri->normals = (float*)arena->alloc(sizeof(Vec3f)*tri->vertices_n);
  tri->indices = (uint32_t*)arena->alloc(3 * sizeof(uint32_t)*tri->triangles_n);

  float h2 = 0.5f*cy.height;
  float h[2] = { -h2, h2 };
  unsigned o = 0;
  auto * V = (Vec3f*)tri->vertices;
  auto * N = (Vec3f*)tri->normals;
  auto * T = (Vec2f*)tri->texCoords;
  auto * I = tri->indices;
  auto * S = tri->smoothingGroups;

  if (texcoords == false) {
    if (shell) {
      for (unsigned i = 0; i < samples; i++) {
        for (unsigned k = 0; k < 2; k++) {
          *V++ = Vec3f(cosSinRadius[i], h[k]);
          *N++ = Vec3f(cosSin[i], 0);
        }
      }
    }
    for (unsigned k = 0; k < 2; k++) {
      if (cap[k] == false) continue;
      for (unsigned i = 0; i < samples; i++) {
        *V++ = Vec3f(cosSinRadius[i], h[k]);
        *N++ = Vec3f(0, 0, -1);
      }
    }
    if (shell) {
      for (unsigned i = 0; i < samples; i++) {
        unsigned ii = (i + 1) % samples;
        I = quadIndices_(I, o, 2 * i, 2 * ii, 2 * ii + 1, 2 * i + 1);
      }
      o += 2 * samples;
      if (smoothingGroups) {
        for (unsigned i = 0; i < 2 * samples; i++) {
          *S++ = geo->componentId;
        }
      }
    }
    for (unsigned k = 0; k < 2; k++) {
      if (cap[k] == false) continue;
      I = tessellateCircle_(I, o, k, samples);
      o += samples;
      if (smoothingGroups) {
        auto sgrp = newSmoothingGroup();
        for (unsigned i = 0; i < (samples - 2); i++) {
          *S++ = sgrp;
        }
      }
    }
  }
  else {
    auto circumference = scale * twopi * cy.radius;
    if (shell) {
      float distances[2] = {
        geo->distances[0],
        geo->distances[1]
      };
      auto iScale = distances[0] < distances[1] ? 1.f / samples : -1.f / samples;
      auto iConst = distances[0] < distances[1] ? 0 : 1;

      for (unsigned i = 0; i < samples; i++) {
        for (unsigned k = 0; k < 2; k++) {
          *V++ = Vec3f(cosSinRadius[i], h[k]);
          *N++ = Vec3f(cosSin[i], 0);
          *T++ = encodeCylindricalTexCoord(iScale * i + iConst, distances[k], circumference);
        }
      }
      for (unsigned k = 0; k < 2; k++) {// clip texture seam.
        *V++ = Vec3f(cosSinRadius[0], h[k]);
        *N++ = Vec3f(cosSin[0], 0);
        *T++ = encodeCylindricalTexCoord(1.f- iConst, distances[k], circumference);
      }
    }
    for (unsigned j = 0; j < 2; j++) {
      if (cap[j] == false) continue;
      for (unsigned i = 0; i < samples; i++) {
        *V++ = Vec3f(cosSinRadius[i], h[j]);
        *N++ = Vec3f(0, 0, j == 0 ? -1.f : 1.f);
        *T++ = encodeCircularTexCoord(cosSin[i], scale * cy.radius);
      }
    }
    if (shell) {
      for (unsigned i = 0; i < samples; i++) {
        unsigned ii = (i + 1);
        I = quadIndices_(I, o, 2 * i, 2 * ii, 2 * ii + 1, 2 * i + 1);
      }
      o += 2 * (samples + 1);
      if (smoothingGroups) {
        for (unsigned i = 0; i < 2 * samples; i++) {
          *S++ = geo->componentId;
        }
      }
    }
    for (unsigned k = 0; k < 2; k++) {
      if (cap[k] == false) continue;
      I = tessellateCircle_(I, o, k, samples);
      o += samples;
      if (smoothingGroups) {
        auto sgrp = newSmoothingGroup();
        for (unsigned i = 0; i < (samples - 2); i++) {
          *S++ = sgrp;
        }
      }
    }
    assert((T - (Vec2f*)tri->texCoords) == tri->vertices_n);
  }

  assert((V - (Vec3f*)tri->vertices) == tri->vertices_n);
  assert((N - (Vec3f*)tri->normals) == tri->vertices_n);
  assert((I -tri->indices) == 3 * tri->triangles_n);
  assert(o == tri->vertices_n);
  if (smoothingGroups) {
    assert((S - tri->smoothingGroups) == tri->triangles_n);
  }
  return tri;
}


Triangulation* TriangulationFactory::sphereBasedShape(Arena* arena, const Geometry* geo, float radius, float arc, float shift_z, float scale_z, float scale)
{
  unsigned segments = sagittaBasedSegmentCount(twopi, radius, scale);
  unsigned samples = segments;  // Assumed to be closed

  auto * tri = arena->alloc<Triangulation>();
  tri->error = sagittaBasedError(twopi, radius, scale, samples);

  bool is_sphere = false;
  if (pi - 1e-3 <= arc) {
    arc = pi;
    is_sphere = true;
  }

  unsigned min_rings = 3;// arc <= half_pi ? 2 : 3;
  unsigned rings = unsigned(std::max(float(min_rings), scale_z * samples*arc*(1.f / twopi)));

  u0.resize(rings);
  t0.resize(2 * rings);
  auto theta_scale = arc / (rings - 1);
  for (unsigned r = 0; r < rings; r++) {
    float theta = theta_scale * r;
    t0[2 * r + 0] = std::cos(theta);
    t0[2 * r + 1] = std::sin(theta);
    u0[r] = unsigned(std::max(3.f, t0[2 * r + 1] * samples));  // samples in this ring
  }
  u0[0] = 1;
  if (is_sphere) {
    u0[rings - 1] = 1;
  }

  unsigned s = 0;
  for (unsigned r = 0; r < rings; r++) {
    s += u0[r];
  }


  tri->vertices_n = s;
  tri->vertices = (float*)arena->alloc(3 * sizeof(float)*tri->vertices_n);
  tri->normals = (float*)arena->alloc(3 * sizeof(float)*tri->vertices_n);

  unsigned l = 0;
  for (unsigned r = 0; r < rings; r++) {
    auto nz = t0[2 * r + 0];
    auto z = radius * scale_z * nz + shift_z;
    auto w = t0[2 * r + 1];
    auto n = u0[r];

    auto phi_scale = twopi / n;
    for (unsigned i = 0; i < n; i++) {
      auto phi = phi_scale * i + geo->sampleStartAngle;
      auto nx = w * std::cos(phi);
      auto ny = w * std::sin(phi);
      l = vertex(tri->normals, tri->vertices, l, nx, ny, nz / scale_z, radius*nx, radius*ny, z);
    }
  }
  assert(l == 3 * tri->vertices_n);

  unsigned o_c = 0;
  indices.clear();
  for (unsigned r = 0; r + 1 < rings; r++) {
    auto n_c = u0[r];
    auto n_n = u0[r + 1];
    auto o_n = o_c + n_c;

    if (n_c < n_n) {
      for (unsigned i_n = 0; i_n < n_n; i_n++) {
        unsigned ii_n = (i_n + 1);
        unsigned i_c = (n_c*(i_n + 1)) / n_n;
        unsigned ii_c = (n_c*(ii_n + 1)) / n_n;

        i_c %= n_c;
        ii_c %= n_c;
        ii_n %= n_n;

        if (i_c != ii_c) {
          indices.push_back(o_c + i_c);
          indices.push_back(o_n + ii_n);
          indices.push_back(o_c + ii_c);
        }
        assert(i_n != ii_n);
        indices.push_back(o_c + i_c);
        indices.push_back(o_n + i_n);
        indices.push_back(o_n + ii_n);
      }
    }
    else {
      for (unsigned i_c = 0; i_c < n_c; i_c++) {
        auto ii_c = (i_c + 1);
        unsigned i_n = (n_n*(i_c + 0)) / n_c;
        unsigned ii_n = (n_n*(ii_c + 0)) / n_c;

        i_n %= n_n;
        ii_n %= n_n;
        ii_c %= n_c;

        assert(i_c != ii_c);
        indices.push_back(o_c + i_c);
        indices.push_back(o_n + ii_n);
        indices.push_back(o_c + ii_c);

        if (i_n != ii_n) {
          indices.push_back(o_c + i_c);
          indices.push_back(o_n + i_n);
          indices.push_back(o_n + ii_n);
        }
      }
    }
    o_c = o_n;
  }

  tri->triangles_n = unsigned(indices.size() / 3);
  tri->indices = (uint32_t*)arena->dup(indices.data(), 3 * sizeof(uint32_t)*tri->triangles_n);

  return tri;
}


Triangulation* TriangulationFactory::facetGroup(Arena* arena, const Geometry* geo, float scale)
{
  auto & fg = geo->facetGroup;

  vertices.clear();
  normals.clear();
  indices.clear();
  for (unsigned p = 0; p < fg.polygons_n; p++) {
    auto & poly = fg.polygons[p];

    if (poly.contours_n == 1 && poly.contours[0].vertices_n == 3) {
      auto & cont = poly.contours[0];
      auto vo = uint32_t(vertices.size()) / 3;

      vertices.resize(vertices.size() + 3 * 3);
      normals.resize(vertices.size());

      std::memcpy(vertices.data() + 3 * vo, cont.vertices, 3 * 3 * sizeof(float));
      std::memcpy(normals.data() + 3 * vo, cont.normals, 3 * 3 * sizeof(float));

      indices.push_back(vo + 0);
      indices.push_back(vo + 1);
      indices.push_back(vo + 2);
    }
    else if (poly.contours_n == 1 && poly.contours[0].vertices_n == 4) {
      auto & cont = poly.contours[0];
      auto & V = cont.vertices;
      auto vo = uint32_t(vertices.size()) / 3;

      vertices.resize(vertices.size() + 3 * 4);
      normals.resize(vertices.size());

      std::memcpy(vertices.data() + 3 * vo, cont.vertices, 3 * 4 * sizeof(float));
      std::memcpy(normals.data() + 3 * vo, cont.normals, 3 * 4 * sizeof(float));

      // find least folding diagonal

      float v01[3], v12[3], v23[3], v30[3];
      sub3(v01, V + 3 * 1, V + 3 * 0);
      sub3(v12, V + 3 * 2, V + 3 * 1);
      sub3(v23, V + 3 * 3, V + 3 * 2);
      sub3(v30, V + 3 * 0, V + 3 * 3);

      float n0[3], n1[3], n2[3], n3[3];
      cross3(n0, v01, v30);
      cross3(n1, v12, v01);
      cross3(n2, v23, v12);
      cross3(n3, v30, v23);

      if (dot3(n0, n2) < dot3(n1, n3)) {
        indices.push_back(vo + 0);
        indices.push_back(vo + 1);
        indices.push_back(vo + 2);

        indices.push_back(vo + 2);
        indices.push_back(vo + 3);
        indices.push_back(vo + 0);
      }
      else {
        indices.push_back(vo + 3);
        indices.push_back(vo + 0);
        indices.push_back(vo + 1);

        indices.push_back(vo + 1);
        indices.push_back(vo + 2);
        indices.push_back(vo + 3);
      }
    }
    else  {

      bool anyData = false;

      BBox3f bbox = createEmptyBBox3f();
      for (unsigned c = 0; c < poly.contours_n; c++) {
        for (unsigned i = 0; i < poly.contours[c].vertices_n; i++) {
          const auto p = Vec3f(poly.contours[c].vertices + 3 * i);
          engulf(bbox, p);
        }
      }
      auto m = 0.5f*(Vec3f(bbox.min) + Vec3f(bbox.max));

      auto tess = tessNewTess(nullptr);
      for (unsigned c = 0; c < poly.contours_n; c++) {
        auto & cont = poly.contours[c];
        if (cont.vertices_n < 3) {
          logger(1, "Ignoring degenerate contour with %d vertices.", cont.vertices_n);
          continue;
        }
        vec3.resize(cont.vertices_n);
        for (unsigned i = 0; i < cont.vertices_n; i++) {
          vec3[i] = Vec3f(cont.vertices + 3 * i) - m;
        }
        tessAddContour(tess, 3, vec3.data(), 3 * sizeof(float), cont.vertices_n);
        //tessAddContour(tess, 3, cont.vertices, 3 * sizeof(float), cont.vertices_n);
        anyData = true;
      }

      if (anyData == false) {
        logger(1, "Ignoring polygon with no valid contours.");
      }
      else {
        if (tessTesselate(tess, TESS_WINDING_ODD, TESS_POLYGONS, 3, 3, nullptr)) {
          auto vo = uint32_t(vertices.size()) / 3;
          auto vn = unsigned(tessGetVertexCount(tess));

          vertices.resize(vertices.size() + 3 * vn);

          auto * src = tessGetVertices(tess);
          for (unsigned i = 0; i < vn; i++) {
            auto p = Vec3f((float*)(src + 3 * i)) + m;
            write(vertices.data() + 3 * (vo + i), p);
          }

          //std::memcpy(vertices.data() + 3 * vo, tessGetVertices(tess), 3 * vn * sizeof(float));

          auto * remap = tessGetVertexIndices(tess);
          normals.resize(vertices.size());
          for (unsigned i = 0; i < vn; i++) {
            if (remap[i] != TESS_UNDEF) {
              unsigned ix = remap[i];
              for (unsigned c = 0; c < poly.contours_n; c++) {
                auto & cont = poly.contours[c];
                if (ix < cont.vertices_n) {
                  normals[3 * (vo + i) + 0] = cont.normals[3 * ix + 0];
                  normals[3 * (vo + i) + 1] = cont.normals[3 * ix + 1];
                  normals[3 * (vo + i) + 2] = cont.normals[3 * ix + 2];
                  break;
                }
                ix -= cont.vertices_n;
              }
            }
          }

          auto io = uint32_t(indices.size());
          auto * elements = tessGetElements(tess);
          auto elements_n = unsigned(tessGetElementCount(tess));
          for (unsigned e = 0; e < elements_n; e++) {
            auto ix = elements + 3 * e;
            if ((ix[0] != TESS_UNDEF) && (ix[1] != TESS_UNDEF) && (ix[2] != TESS_UNDEF)) {
              indices.push_back(ix[0] + vo);
              indices.push_back(ix[1] + vo);
              indices.push_back(ix[2] + vo);
            }
          }
        }
      }

      tessDeleteTess(tess);
    }
  }

  assert(vertices.size() == normals.size());

  Triangulation* tri = arena->alloc<Triangulation>();
  tri->error = 0.f;

  if (!indices.empty()) {
    tri->vertices_n = uint32_t(vertices.size() / 3);
    tri->triangles_n =  uint32_t(indices.size() / 3);

    tri->vertices = (float*)arena->dup(vertices.data(), sizeof(float) * 3 * tri->vertices_n);
    tri->normals = (float*)arena->dup(normals.data(), sizeof(float) * 3 * tri->vertices_n);
    tri->indices = (uint32_t*)arena->dup(indices.data(), sizeof(uint32_t) * 3 * tri->triangles_n);
  }

  return tri;
}

