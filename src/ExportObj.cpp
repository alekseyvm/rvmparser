#include <cassert>
#include <cstdio>
#include <string>
#include <algorithm>
#include "ExportObj.h"
#include "Store.h"
#include "LinAlgOps.h"

namespace {

  bool open_w(FILE** f, const char* path)
  {
    auto err = fopen_s(f, path, "w");
    if (err == 0) return true;

    char buf[1024];
    if (strerror_s(buf, sizeof(buf), err) != 0) {
      buf[0] = '\0';
    }
    fprintf(stderr, "Failed to open %s for writing: %s", path, buf);
    return false;
  }

  void wireBoundingBox(FILE* out, unsigned& off_v, const BBox3f& bbox)
  {
    for (unsigned i = 0; i < 8; i++) {
      float px = (i & 1) ? bbox.min[0] : bbox.min[3];
      float py = (i & 2) ? bbox.min[1] : bbox.min[4];
      float pz = (i & 4) ? bbox.min[2] : bbox.min[5];
      fprintf(out, "v %f %f %f\n", px, py, pz);
    }
    fprintf(out, "l %d %d %d %d %d\n",
            off_v + 0, off_v + 1, off_v + 3, off_v + 2, off_v + 0);
    fprintf(out, "l %d %d %d %d %d\n",
            off_v + 4, off_v + 5, off_v + 7, off_v + 6, off_v + 4);
    fprintf(out, "l %d %d\n", off_v + 0, off_v + 4);
    fprintf(out, "l %d %d\n", off_v + 1, off_v + 5);
    fprintf(out, "l %d %d\n", off_v + 2, off_v + 6);
    fprintf(out, "l %d %d\n", off_v + 3, off_v + 7);
    off_v += 8;
  }

}


ExportObj::~ExportObj()
{
  if (out) {
    fclose(out);
  }
  if (mtl) {
    fclose(mtl);
  }
}

bool ExportObj::open(const char* path_obj, const char* path_mtl)
{
  if (!open_w(&out, path_obj)) return false;
  if (!open_w(&mtl, path_mtl)) return false;

  std::string mtllib(path_mtl);
  auto l = mtllib.find_last_of("/\\");
  if (l != std::string::npos) {
    mtllib = mtllib.substr(l + 1);
  }

  fprintf(out, "mtllib %s\n", mtllib.c_str());

  if (groupBoundingBoxes) {
    fprintf(mtl, "newmtl group_bbox\n");
    fprintf(mtl, "Ka 1 0 0\n");
    fprintf(mtl, "Kd 1 0 0\n");
    fprintf(mtl, "Ks 0 0 0\n");
  }

  return true;
}


void ExportObj::init(class Store& store)
{
  assert(out);
  assert(mtl);

  this->store = &store;
  conn = store.conn;

  stack.accommodate(store.groupCountAllocated());
  for (auto * line = store.getFirstDebugLine(); line != nullptr; line = line->next) {
    useColor(nullptr, line->color);
    fprintf(out, "v %f %f %f\n", line->a[0], line->a[1], line->a[2]);
    fprintf(out, "v %f %f %f\n", line->b[0], line->b[1], line->b[2]);
    fprintf(out, "l -1 -2\n");
    off_v += 2;
  }
}

void ExportObj::beginFile(Group* group)
{
  fprintf(out, "# %s\n", group->file.info);
  fprintf(out, "# %s\n", group->file.note);
  fprintf(out, "# %s\n", group->file.date);
  fprintf(out, "# %s\n", group->file.user);
}

void ExportObj::endFile() {

  //fprintf(out, "usemtl red_line\n");

  //auto l = 0.05f;
  //if (anchors && conn) {
  //  for (unsigned i = 0; i < conn->anchor_n; i++) {
  //    auto * p = conn->p + 3 * i;
  //    auto * n = conn->anchors[i].n;

  //    fprintf(out, "v %f %f %f\n", p[0], p[1], p[2]);
  //    fprintf(out, "v %f %f %f\n", p[0] + l * n[0], p[1] + l * n[1], p[2] + l * n[2]);
  //    fprintf(out, "l %d %d\n", (int)off_v, (int)off_v + 1);
  //    off_v += 2;
  //  }
  //}

  fprintf(out, "# End of file\n");
}

void ExportObj::beginModel(Group* group)
{
  fprintf(out, "# Model project=%s, name=%s\n", group->model.project, group->model.name);
}

void ExportObj::endModel() { }

void ExportObj::beginGroup(Group* group)
{
  for (unsigned i = 0; i < 3; i++) curr_translation[i] = group->group.translation[i];

  stack[stack_p++] = group->group.name;

  fprintf(out, "o %s", stack[0]);
  for (unsigned i = 1; i < stack_p; i++) {
    fprintf(out, "/%s", stack[i]);
  }
  fprintf(out, "\n");
 

//  fprintf(out, "o %s\n", group->group.name);


  if (groupBoundingBoxes && !isEmpty(group->group.bboxWorld)) {
    fprintf(out, "usemtl group_bbox\n");
    wireBoundingBox(out, off_v, group->group.bboxWorld);
  }

}

void ExportObj::EndGroup() {
  assert(stack_p);
  stack_p--;
}

namespace {

  void getMidpoint(Vec3f& p, Geometry* geo)
  {
    switch (geo->kind) {
    case Geometry::Kind::CircularTorus: {
      auto & ct = geo->circularTorus;
      auto c = std::cos(0.5f * ct.angle);
      auto s = std::sin(0.5f * ct.angle);
      p.x = ct.offset * c;
      p.y = ct.offset * s;
      p.z = 0.f;
      break;
    }
    default:
      p = Vec3f(0.f);
      break;
    }
    p = mul(geo->M_3x4, p);
  }

}

void ExportObj::useColor(const char* colorName, uint32_t color)
{
  if (currentColor == color) return;
  currentColor = color;

  char name[9];
  name[0] = '0';
  name[1] = 'x';
  for (unsigned k = 0; k < 6; k++) {
    auto v = (color >> (4 * (5 - k))) & 0xf;
    if (v < 10) name[k+2] = '0' + v;
    else name[k+2] = 'a' + v - 10;
  }
  name[8] = '\0';
  auto * hexName = store->strings.intern(name);

  if (!definedColors.get(uint64_t(hexName))) {
    definedColors.insert(uint64_t(hexName), 1);

    auto r = (1.f / 255.f)*((color >> 16) & 0xFF);
    auto g = (1.f / 255.f)*((color >> 8) & 0xFF);
    auto b = (1.f / 255.f)*((color) & 0xFF);


    fprintf(mtl, "newmtl %s\n", hexName);
    fprintf(mtl, "Ka %f %f %f\n", (2.f / 3.f)*r, (2.f / 3.f)*g, (2.f / 3.f)*b);
    fprintf(mtl, "Kd %f %f %f\n", r, g, b);
    fprintf(mtl, "Ks 0.5 0.5 0.5\n");
  }
  fprintf(out, "usemtl %s\n", hexName);

}

void ExportObj::geometry(struct Geometry* geometry)
{
  

  const auto & M = geometry->M_3x4;


  //if (geometry->colorName == nullptr) {
  //  geometry->colorName = store->strings.intern("default");
  //}

  //if (geometry->kind == Geometry::Kind::Box) {
  //  geometry->colorName = store->strings.intern("blah-red");
  //  geometry->color = 0x880088;
  //}
  //if (geometry->kind == Geometry::Kind::Pyramid) {
  //  geometry->colorName = store->strings.intern("blah-green");
  //  geometry->color = 0x008800;
  //}
  //if (geometry->kind == Geometry::Kind::RectangularTorus) {
  //  geometry->colorName = store->strings.intern("blah-blue");
  //  geometry->color = 0x000088;
  //}
  //if (geometry->kind == Geometry::Kind::FacetGroup) {
  //  geometry->colorName = store->strings.intern("blah-redgg");
  //  geometry->color = 0x888800;
  //}


  

  auto scale = 1.f;
  
  if (geometry->kind == Geometry::Kind::Line) {
    useColor(geometry->colorName, geometry->color);

    auto a = scale * mul(geometry->M_3x4, Vec3f(geometry->line.a, 0, 0));
    auto b = scale * mul(geometry->M_3x4, Vec3f(geometry->line.b, 0, 0));
    fprintf(out, "v %f %f %f\n", a.x, a.y, a.z);
    fprintf(out, "v %f %f %f\n", b.x, b.y, b.z);
    fprintf(out, "l -1 -2\n");
    off_v += 2;
  }
  else {
    assert(geometry->triangulation);
    auto * tri = geometry->triangulation;

    if (tri->dPdu && tri->dPdv) {

      for (unsigned i = 0; i < tri->vertices_n; i++) {
        auto & uu = tri->dPdu[i];

        if (uu.x == 0.f && uu.y == 0.f && uu.z == 0.f) continue;
        continue;

        auto p = scale * mul(geometry->M_3x4, Vec3f(tri->vertices + 3 * i));
        auto n = normalize(mul(Mat3f(geometry->M_3x4.data), Vec3f(tri->normals + 3 * i)));
        auto u = normalize(mul(Mat3f(geometry->M_3x4.data), tri->dPdu[i]));
        auto v = normalize(mul(Mat3f(geometry->M_3x4.data), tri->dPdv[i]));

        auto o = p + 0.0001f*n;
        auto ou = o + 0.05f*u;
        auto ov = o + 0.05f*v;
        auto on = o + 0.05f*cross(u, v);

        fprintf(out, "v %f %f %f\n", o.x, o.y, o.z);
        fprintf(out, "v %f %f %f\n", ou.x, ou.y, ou.z);
        fprintf(out, "v %f %f %f\n", ov.x, ov.y, ov.z);
        fprintf(out, "v %f %f %f\n", on.x, on.y, on.z);

        useColor(nullptr, 0xff0000u);
        fprintf(out, "l -3 -4\n");
        useColor(nullptr, 0x00ff00u);
        fprintf(out, "l -2 -4\n");
        useColor(nullptr, 0x0000ffu);
        fprintf(out, "l -1 -4\n");
        off_v += 4;
      }
    }

    useColor(geometry->colorName, geometry->color);

    if (tri->indices != 0) {
      //fprintf(out, "g\n");
      if (geometry->triangulation->error != 0.f) {
        fprintf(out, "# error=%f\n", geometry->triangulation->error);
      }
      for (size_t i = 0; i < 3 * tri->vertices_n; i += 3) {

        auto p = scale * mul(geometry->M_3x4, Vec3f(tri->vertices + i));
        auto n = normalize(mul(Mat3f(geometry->M_3x4.data), Vec3f(tri->normals + i)));

        fprintf(out, "v %f %f %f\n", p.x, p.y, p.z);
        fprintf(out, "vn %f %f %f\n", n.x, n.y, n.z);
      }

      if (tri->texCoords) {
        for (size_t i = 0; i < tri->vertices_n; i++) {
          const Vec2f vt(tri->texCoords + 2 * i);
          fprintf(out, "vt %f %f\n", vt.x, vt.y);
        }
      }
      else {
        for (size_t i = 0; i < tri->vertices_n; i++) {
          auto p = scale * mul(geometry->M_3x4, Vec3f(tri->vertices + 3 * i));
          fprintf(out, "vt %f %f\n", 0 * p.x, 0 * p.y);
        }
      }

      uint32_t currentSmoothingGroup = 0;
      for (size_t i = 0; i < tri->triangles_n; i++) {
        auto smoothingGroup = tri->smoothingGroups ? tri->smoothingGroups[i] : 0;
        if (i == 0 || currentSmoothingGroup != smoothingGroup) {
          currentSmoothingGroup = smoothingGroup;
          if (smoothingGroup == 0) {
            fprintf(out, "s off\n");
          }
          else {
            fprintf(out, "s %u\n", smoothingGroup);
          }
        }
        auto a = tri->indices[3 * i + 0];
        auto b = tri->indices[3 * i + 1];
        auto c = tri->indices[3 * i + 2];
        fprintf(out, "f %d/%d/%d %d/%d/%d %d/%d/%d\n",
                a + off_v, a + off_t, a + off_n,
                b + off_v, b + off_t, b + off_n,
                c + off_v, c + off_t, c + off_n);
      }
    
      off_v += tri->vertices_n;
      off_n += tri->vertices_n;
      off_t += tri->vertices_n;
    }
  }

  //if (primitiveBoundingBoxes) {
  //  fprintf(out, "usemtl magenta\n");

  //  for (unsigned i = 0; i < 8; i++) {
  //    float px = (i & 1) ? geometry->bbox[0] : geometry->bbox[3];
  //    float py = (i & 2) ? geometry->bbox[1] : geometry->bbox[4];
  //    float pz = (i & 4) ? geometry->bbox[2] : geometry->bbox[5];

  //    float Px = M[0] * px + M[3] * py + M[6] * pz + M[9];
  //    float Py = M[1] * px + M[4] * py + M[7] * pz + M[10];
  //    float Pz = M[2] * px + M[5] * py + M[8] * pz + M[11];

  //    fprintf(out, "v %f %f %f\n", Px, Py, Pz);
  //  }
  //  fprintf(out, "l %d %d %d %d %d\n",
  //          off_v + 0, off_v + 1, off_v + 3, off_v + 2, off_v + 0);
  //  fprintf(out, "l %d %d %d %d %d\n",
  //          off_v + 4, off_v + 5, off_v + 7, off_v + 6, off_v + 4);
  //  fprintf(out, "l %d %d\n", off_v + 0, off_v + 4);
  //  fprintf(out, "l %d %d\n", off_v + 1, off_v + 5);
  //  fprintf(out, "l %d %d\n", off_v + 2, off_v + 6);
  //  fprintf(out, "l %d %d\n", off_v + 3, off_v + 7);
  //  off_v += 8;
  //}


  //for (unsigned k = 0; k < 6; k++) {
  //  auto other = geometry->conn_geo[k];
  //  if (geometry < other) {
  //    fprintf(out, "usemtl blue_line\n");
  //    float p[3];
  //    getMidpoint(p, geometry);
  //    fprintf(out, "v %f %f %f\n", p[0], p[1], p[2]);
  //    getMidpoint(p, other);
  //    fprintf(out, "v %f %f %f\n", p[0], p[1], p[2]);
  //    fprintf(out, "l %d %d\n", off_v, off_v + 1);

  //    off_v += 2;
  //  }
  //}

}
