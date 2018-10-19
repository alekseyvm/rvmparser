#pragma once
#include "Common.h"
#include "Store.h"

struct Interface
{
  enum struct Kind {
    Undefined,
    Square,
    Circular
  };
  Kind kind = Kind::Undefined;

  union {
    struct {
      Vec3f p[4];
    } square;
    struct {
      float radius;
    } circular;
  };

};

Interface getInterface(const Geometry* geo, unsigned o);

bool doInterfacesMatch(const Geometry* geo, const Connection* con);
