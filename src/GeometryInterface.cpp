#include <cassert>
#include "GeometryInterface.h"
#include "LinAlgOps.h"

Interface getInterface(const Geometry* geo, unsigned o)
{
  Interface interface;
  auto * connection = geo->connections[o];
  auto ix = connection->geo[0] == geo ? 0 : 1;
  auto scale = getScale(geo->M_3x4);
  switch (geo->kind) {
  case Geometry::Kind::Pyramid: {
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
        Vec3f(-bx - ox, -by - oy, -h2),
        Vec3f(bx - ox, -by - oy, -h2),
        Vec3f(bx - ox,  by - oy, -h2),
        Vec3f(-bx - ox,  by - oy, -h2)
      },
      {
         Vec3f(-tx + ox, -ty + oy, h2),
         Vec3f(tx + ox, -ty + oy, h2),
         Vec3f(tx + ox,  ty + oy, h2),
         Vec3f(-tx + ox,  ty + oy, h2)
      },
    };

    interface.kind = Interface::Kind::Square;
    if (o < 4) {
      unsigned oo = (o + 1) & 3;
      interface.square.p[0] = mul(geo->M_3x4, quad[0][o]);
      interface.square.p[1] = mul(geo->M_3x4, quad[0][oo]);
      interface.square.p[2] = mul(geo->M_3x4, quad[1][oo]);
      interface.square.p[3] = mul(geo->M_3x4, quad[1][o]);
    }
    else {
      for (unsigned k = 0; k < 4; k++) interface.square.p[k] = mul(geo->M_3x4, quad[o - 4][k].data);
    }
    break;
  }
  case Geometry::Kind::Box: {
    auto & box = geo->box;
    auto xp = 0.5f * box.lengths[0]; auto xm = -xp;
    auto yp = 0.5f * box.lengths[1]; auto ym = -yp;
    auto zp = 0.5f * box.lengths[2]; auto zm = -zp;
    Vec3f V[6][4] = {
      { Vec3f(xm, ym, zp), Vec3f(xm, yp, zp), Vec3f(xm, yp, zm), Vec3f(xm, ym, zm) },
      { Vec3f(xp, ym, zm), Vec3f(xp, yp, zm), Vec3f(xp, yp, zp), Vec3f(xp, ym, zp) },
      { Vec3f(xp, ym, zm), Vec3f(xp, ym, zp), Vec3f(xm, ym, zp), Vec3f(xm, ym, zm) },
      { Vec3f(xm, yp, zm), Vec3f(xm, yp, zp), Vec3f(xp, yp, zp), Vec3f(xp, yp, zm) },
      { Vec3f(xm, yp, zm), Vec3f(xp, yp, zm), Vec3f(xp, ym, zm), Vec3f(xm, ym, zm) },
      { Vec3f(xm, ym, zp), Vec3f(xp, ym, zp), Vec3f(xp, yp, zp), Vec3f(xm, yp, zp) }
    };
    for (unsigned k = 0; k < 4; k++) interface.square.p[k] = mul(geo->M_3x4, V[o][k].data);
    break;
  }
  case Geometry::Kind::RectangularTorus: {
    auto & tor = geo->rectangularTorus;
    auto h2 = 0.5f*tor.height;
    float square[4][2] = {
      { tor.outer_radius, -h2 },
      { tor.inner_radius, -h2 },
      { tor.inner_radius, h2 },
      { tor.outer_radius, h2 },
    };
    if (o == 0) {
      for (unsigned k = 0; k < 4; k++) {
        interface.square.p[k] = mul(geo->M_3x4, Vec3f(square[k][0], 0.f, square[k][1]));
      }
    }
    else {
      for (unsigned k = 0; k < 4; k++) {
        interface.square.p[k] = mul(geo->M_3x4, Vec3f(square[k][0] * cos(tor.angle),
                                                      square[k][0] * sin(tor.angle),
                                                      square[k][1]));
      }
    }
    break;
  }
  case Geometry::Kind::CircularTorus:
    interface.kind = Interface::Kind::Circular;
    interface.circular.radius = scale * geo->circularTorus.radius;
    break;

  case Geometry::Kind::EllipticalDish:
    interface.kind = Interface::Kind::Circular;
    interface.circular.radius = scale * 0.5f*geo->ellipticalDish.diameter;
    break;

  case Geometry::Kind::SphericalDish: {
    float r_circ = 0.5f * geo->sphericalDish.diameter;
    auto h = geo->sphericalDish.height;
    float r_sphere = (r_circ*r_circ + h * h) / (2.f*h);
    interface.kind = Interface::Kind::Circular;
    interface.circular.radius = scale * r_sphere;
    break;
  }
  case Geometry::Kind::Snout:
    interface.kind = Interface::Kind::Circular;
    interface.circular.radius = scale * (connection->offset[ix] == 0 ? geo->snout.radius_b : geo->snout.radius_t);
    break;
  case Geometry::Kind::Cylinder:
    interface.kind = Interface::Kind::Circular;
    interface.circular.radius = scale * geo->cylinder.radius;
    break;
  case Geometry::Kind::Sphere:
  case Geometry::Kind::Line:
  case Geometry::Kind::FacetGroup:
    interface.kind = Interface::Kind::Undefined;
    break;

  default:
    assert(false && "Unhandled primitive type");
    break;
  }
  return interface;
}

bool doInterfacesMatch(const Geometry* geo, const Connection* con, float radialEpsilon)
{
  bool isFirst = geo == con->geo[0];

  auto * thisGeo = con->geo[isFirst ? 0 : 1];
  auto thisOffset = con->offset[isFirst ? 0 : 1];
  auto thisIFace = getInterface(thisGeo, thisOffset);

  auto * thatGeo = con->geo[isFirst ? 1 : 0];
  auto thatOffset = con->offset[isFirst ? 1 : 0];
  auto thatIFace = getInterface(thatGeo, thatOffset);


  if (thisIFace.kind != thatIFace.kind) return false;

  if (thisIFace.kind == Interface::Kind::Circular) {

    return std::abs(thisIFace.circular.radius - thatIFace.circular.radius) < radialEpsilon;

  }
  else {
    for (unsigned j = 0; j < 4; j++) {
      bool found = false;
      for (unsigned i = 0; i < 4; i++) {
        if (distanceSquared(thisIFace.square.p[j], thatIFace.square.p[i]) < radialEpsilon*radialEpsilon) {
          found = true;
        }
      }
      if (!found) return false;
    }
    return true;
  }
}
