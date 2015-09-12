#include "common.cl"

typedef struct {
  float3 position, next_velocity;
  int collision_happened;
  float time_elapsed;
} collision_response;

typedef struct {
  float3 collision_point, surface_normal;
  float penetration_depth;
  int collision_happened;
} collision;

int detect_collision(collision *c, float3 p0, float3 p1,
                     global const float *face_normals,
                     global const float *vertices, global const uint *indices,
                     uint face_count) {
  c->collision_happened = 0;

  for (uint i = 0; i < face_count; ++i) {
    float3 normal = {
        face_normals[3 * i + 0], face_normals[3 * i + 1],
        face_normals[3 * i + 2],
    };

    if (dot(normal, p1 - p0) / (length(normal) * length(p1 - p0)) <= 0) {
      normal = -normal;
    }

    float3 v0 = {
        vertices[3 * indices[3 * i + 0] + 0],
        vertices[3 * indices[3 * i + 0] + 1],
        vertices[3 * indices[3 * i + 0] + 2],
    };

    float3 v1 = {
        vertices[3 * indices[3 * i + 1] + 0],
        vertices[3 * indices[3 * i + 1] + 1],
        vertices[3 * indices[3 * i + 1] + 2],
    };

    float3 v2 = {
        vertices[3 * indices[3 * i + 2] + 0],
        vertices[3 * indices[3 * i + 2] + 1],
        vertices[3 * indices[3 * i + 2] + 2],
    };

    float3 u = v1 - v0;
    float3 v = v2 - v0;

    float denom = dot(normal, p1 - p0);

    if (denom == 0.f) {
      continue;
    }

    float r = dot(normal, v0 - p0) / denom;

    if (0 <= r && r <= 1) {
      float3 intersect = p0 + r * (p1 - p0);
      float3 w = intersect - v0;
      float uv, wv, vv, wu, uu;

      uv = dot(u, v);
      wv = dot(w, v);
      vv = dot(v, v);
      wu = dot(w, u);
      uu = dot(u, u);

      float denom = uv * uv - uu * vv;
      float s = (uv * wv - vv * wu) / denom;
      float t = (uv * wu - uu * wv) / denom;

      // Collision
      if (s >= 0 && t >= 0 && s + t <= 1) {
        if (c->collision_happened &&
            length(p0 - intersect) > length(p0 - c->collision_point)) {
          continue;
        }
        c->surface_normal = normal;
        c->collision_point = intersect;
        c->penetration_depth = length(p1 - intersect);
        c->collision_happened = 1;
      }
    }
  }
  return c->collision_happened;
}

int respond(collision_response *response, collision c, float restitution,
            float time_elapsed) {
  // hack to avoid points directly on the faces, the collision detection code
  // should be
  response->position = c.collision_point - (c.surface_normal * 0.001);

  response->next_velocity -=
      (1.f +
       restitution * c.penetration_depth /
           (time_elapsed * length(response->next_velocity))) *
      dot(response->next_velocity, c.surface_normal) * c.surface_normal;

  return 1;
}

// After the collision response, the particle's position is only partially
// updated,
// the advection  must be recursive until the entire movement is completed
collision_response handle_collisions(
    float3 old_position, float3 position, float3 next, float restitution,
    float time_elapsed, global const float *face_normals,
    global const float *vertices, global const uint *indices, uint face_count) {
  collision_response response = {
      position, next, 0, time_elapsed,
  };

  collision c;

  if (detect_collision(&c, old_position, position, face_normals, vertices,
                       indices, face_count)) {
    response.collision_happened = 1;
    respond(&response, c, restitution, time_elapsed);
    response.time_elapsed =
        time_elapsed * (length(response.position - old_position) /
                        length(position - old_position));
  }

  return response;
}