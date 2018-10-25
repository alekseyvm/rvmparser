#include <cassert>
#include <chrono>
#include <algorithm>
#include "Common.h"
#include "Store.h"
#include "GeometryInterface.h"
#include "LinAlgOps.h"

namespace {

  struct StackItem
  {
    Geometry* geo;
    Connection* from;
    float distance;
  };

  struct Context
  {
    Logger logger = nullptr;
    Store* store = nullptr;

    Buffer<Geometry*> geos;
    unsigned geosCount = 0;

    Buffer<StackItem> stack;

  };

  void extractGeometries(Context* context, Group* group)
  {
    for (auto * child = group->groups.first; child; child = child->next) {
      extractGeometries(context, child);
    }
    if (group->kind == Group::Kind::Group) {
      for (auto * geo = group->group.geometries.first; geo; geo = geo->next) {
        context->geos[context->geosCount++] = geo;
      }
    }
  }

  void growFromSeedGeometry(Context* context, Connection* from, Geometry* seed, uint32_t componentId, float distance)
  {
    auto & stack = context->stack;

    if (from) from->temp = componentId;
    assert(seed->componentId == 0);

    uint32_t stackBack = 0;
    uint32_t stackFront = 0;
    stack[stackBack++] = { seed, from, distance };
    while (stackFront < stackBack) {
      auto item = stack[stackFront++];
      auto * thisGeo = item.geo;

      if (thisGeo->componentId != 0) {
        assert(thisGeo->componentId == componentId);
        continue;
      }
      thisGeo->componentId = componentId;

      for (unsigned i = 0; i < 6; i++) {

        auto primitiveLength = 0.f;
        switch (item.geo->kind)
        {
        case Geometry::Kind::CircularTorus:
          primitiveLength = item.geo->circularTorus.offset * item.geo->circularTorus.angle;
          break;
        case Geometry::Kind::Cylinder:
          primitiveLength = item.geo->cylinder.height;
          break;
        default:
          break;
        }
        auto distanceNew = item.distance + getScale(item.geo->M_3x4) * primitiveLength;
        thisGeo->distances[i] = distanceNew;

        auto * con = thisGeo->connections[i];
        if (con == item.from) {
          thisGeo->distances[i] = item.distance;
        }
        else if (con && con->temp == 0) {
          bool isFirst = thisGeo == con->geo[0];
          auto * thatGeo = con->geo[isFirst ? 1 : 0];

          assert(thatGeo != thisGeo);
          assert(con->temp == 0);
          con->temp = componentId;
          stack[stackBack++] = { thatGeo , con, distanceNew };
        }
      }
    }
  }

}


void growComponents(Store* store, Logger logger)
{
  Context context;
  context.logger = logger;
  context.store = store;
  auto time0 = std::chrono::high_resolution_clock::now();
  for (auto * connection = store->getFirstConnection(); connection != nullptr; connection = connection->next) {
    connection->temp = 0;

    for (auto s = 0; s < 2; s++) {
      auto * geo = connection->geo[s];
      unsigned o = ~0u;
      for (unsigned i = 0; i < 6; i++) {
        if (geo->connections[i] == connection) {
          o = i;
        }
      }
      assert(o != ~0u);
      assert(connection->offset[s ? 0 : 1] != ~0u);
    }

  }

  context.stack.accommodate(store->geometryCountAllocated());
  context.geos.accommodate(store->geometryCountAllocated());
  for (auto * root = store->getFirstRoot(); root; root = root->next) {
    extractGeometries(&context, root);
  }

  // First grow open components by picking an end and growing from there.
  unsigned openComponents = 0;
  for (unsigned g = 0; g < context.geosCount; g++) {
    auto * geo = context.geos[g];

    auto unprocessedConnections = 0;
    for (unsigned i = 0; i < 6; i++) {
      if (geo->connections[i] && geo->connections[i]->temp == 0) unprocessedConnections++;
    }

    switch (geo->kind) {
    case Geometry::Kind::Cylinder:
      if (unprocessedConnections == 1) {
        growFromSeedGeometry(&context, nullptr, geo, store->newComponent() , 0.f);
        openComponents++;
      }
      break;
    default:
      break;
    }
  }

  // Then grow closed components by picking arbitrary unprocessed geometries and grow from there.
  unsigned closedComponents = 0;
  for (unsigned g = 0; g < context.geosCount; g++) {
    auto * geo = context.geos[g];
    if (geo->componentId == 0) {
      growFromSeedGeometry(&context, geo->connections[0], geo, store->newComponent(), 0.f);
      closedComponents++;
    }
  }

  auto time1 = std::chrono::high_resolution_clock::now();
  auto e0 = std::chrono::duration_cast<std::chrono::milliseconds>((time1 - time0)).count();

  logger(0, "Processed %d connections -> %d open and %d closed components (%lldms).", store->connectionCountAllocated(), openComponents, closedComponents, e0);
}
