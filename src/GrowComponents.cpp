#include <cassert>
#include <chrono>
#include <algorithm>
#include "Common.h"
#include "Store.h"
#include "GeometryInterface.h"

namespace {

  struct QueueItem
  {
    Geometry* from;
    Connection* connection;
  };

  struct Context
  {
    Logger logger = nullptr;
    Store* store = nullptr;

    Buffer<QueueItem> queue;
    unsigned queueFront = 0;
    unsigned queueBack = 0;

    Buffer<Geometry*> geos;
    unsigned geosCount = 0;

    Buffer<Geometry*> stack;

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

  void growFromSeedGeometry(Context* context, Geometry* seed, uint32_t componentId, float distance)
  {
    auto & stack = context->stack;

    uint32_t stack_p = 0;
    stack[stack_p++] = seed;
    while (stack_p) {
      auto * thisGeo = stack[--stack_p];
      thisGeo->componentId = componentId;
      for (auto * con : thisGeo->connections) {
        if (con == nullptr || con->temp) continue;

        bool isFirst = thisGeo == con->geo[0];
        auto thisOffset = con->offset[isFirst ? 0 : 1];
        auto thisIFace = getInterface(thisGeo, thisOffset);
        auto * thatGeo = con->geo[isFirst ? 1 : 0];
        auto thatOffset = con->offset[isFirst ? 1 : 0];
        auto thatIFace = getInterface(thatGeo, thatOffset);

        con->temp = 1;
        stack[stack_p++] = thatGeo;
      }
    }
  }

}


void growComponents(Store* store, Logger logger)
{
  Context context;
  context.logger = logger;
  context.store = store;
  context.queue.accommodate(store->connectionCountAllocated());
  auto time0 = std::chrono::high_resolution_clock::now();
  for (auto * connection = store->getFirstConnection(); connection != nullptr; connection = connection->next) {
    connection->temp = 0;
  }

  // start at the ends
  context.stack.accommodate(store->geometryCountAllocated());
  context.geos.accommodate(store->geometryCountAllocated());
  for (auto * root = store->getFirstRoot(); root; root = root->next) {
    extractGeometries(&context, root);
  }

  unsigned seedGeometries = 0;
  for (unsigned g = 0; g < context.geosCount; g++) {
    auto * geo = context.geos[g];
    auto unprocessedConnections = 0;
    for (unsigned i = 0; i < 6; i++) {
      if (geo->connections[i] && geo->connections[i]->temp == 0) unprocessedConnections++;
    }

    switch (geo->kind) {
    case Geometry::Kind::Cylinder:
      if (unprocessedConnections == 1) {
        growFromSeedGeometry(&context, geo, store->newComponent() , 0.f);
        seedGeometries++;
      }
      break;
    default:
      break;
    }
  }


  for (auto * connection = store->getFirstConnection(); connection != nullptr; connection = connection->next) {



  }
  auto time1 = std::chrono::high_resolution_clock::now();
  auto e0 = std::chrono::duration_cast<std::chrono::milliseconds>((time1 - time0)).count();

  logger(0, "Processed %d connections %d seed geometries (%lldms).", store->connectionCountAllocated(),seedGeometries, e0);
}
