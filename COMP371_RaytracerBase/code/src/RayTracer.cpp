#include "Sphere.h"
#include "Rectangle.h"
#include "RayTracer.h"
#include "AreaLight.h"
#include "PointLight.h"
#include <random>

using json = nlohmann::json;
using Vector3f = Eigen::Vector3f;
using namespace std;

RayTracer::RayTracer(const json& scene) {
    createGeometry(scene);
    createLights(scene);
    createOutputs(scene);
}


void RayTracer::createGeometry(const json& scene) {
    for (auto& object : scene["geometry"]) {
      if (object["type"] == "sphere") {
        cout<<"Creating sphere"<<endl;

        Sphere sphere(object["centre"][0], object["centre"][1], object["centre"][2], object["radius"]);

        sphere.SetReflectionCoefficients(object["ka"], object["kd"], object["ks"]);

        sphere.SetReflectionColor(object["ac"][0], object["ac"][1], object["ac"][2],
        object["dc"][0], object["dc"][1], object["dc"][2], 
        object["sc"][0], object["sc"][1], object["sc"][2]);

        sphere.SetPhongCoefficients(object["pc"]);

        spheres.push_back(sphere);      
      }
      else if (object["type"] == "rectangle") {
        cout<<"Creating rectangle"<<endl;

          Rectangle rectangle(object["p1"][0], object["p1"][1], object["p1"][2], 
          object["p2"][0], object["p2"][1], object["p2"][2], 
          object["p3"][0], object["p3"][1], object["p3"][2], 
          object["p4"][0], object["p4"][1], object["p4"][2]);

          rectangle.SetReflectionCoefficients(object["ka"], object["kd"], object["ks"]);

          rectangle.SetReflectionColor(object["ac"][0], object["ac"][1], object["ac"][2], 
          object["dc"][0], object["dc"][1], object["dc"][2], 
          object["sc"][0], object["sc"][1], object["sc"][2]);

          rectangle.SetPhongCoefficients(object["pc"]);

          rectangles.push_back(rectangle);
      }
      else {
          cout<<"Error: Unknown object type"<<endl;
      }
    }
}

void RayTracer::createLights(const json& scene) {
    for (auto& light : scene["light"]) {
        if (light["type"] == "point") {
            cout<<"Creating point light"<<endl;

            PointLight pointLight(light["centre"][0], light["centre"][1], light["centre"][2]);

            pointLight.SetLightIntensity(light["id"][0], light["id"][1], light["id"][2], light["is"][0], light["is"][1], light["is"][2]);

            pointLights.push_back(pointLight);
        }
        else if (light["type"] == "area") {
            cout<<"Creating area light"<<endl;

            AreaLight areaLight(light["p1"][0], light["p1"][1], light["p1"][2], 
            light["p2"][0], light["p2"][1], light["p2"][2], 
            light["p3"][0], light["p3"][1], light["p3"][2], 
            light["p4"][0], light["p4"][1], light["p4"][2], light["usecenter"]);

            areaLight.SetLightIntensity(light["id"][0], light["id"][1], light["id"][2], light["is"][0], light["is"][1], light["is"][2]);

            areaLights.push_back(areaLight);
        }
        else {
            cout<<"Error: Unknown light type"<<endl;
        }
    }
}

void RayTracer::createOutputs(const json& scene){
  for (auto& output : scene["output"]) {
    filename = output["filename"];
    size[0] = output["size"][0];
    size[1] = output["size"][1];

    lookat = Vector3f(output["lookat"][0], output["lookat"][1], output["lookat"][2]);
    up = Vector3f(output["up"][0], output["up"][1], output["up"][2]);
    center = Vector3f(output["centre"][0], output["centre"][1], output["centre"][2]);
    bkc = Vector3f(output["bkc"][0], output["bkc"][1], output["bkc"][2]);
    fov = output["fov"];
    try{
      if(output.contains("raysperpixel")){
        raysperpixel[0] = output["raysperpixel"][0];
        raysperpixel[1] = output["raysperpixel"][1];
      }
      if(output.contains("globalillum")){
        globalIllumination = output["globalillum"];
      }
      if(output.contains("maxbounces")){
        maxBounces = output["maxbounces"];
      }
      if(output.contains("probterminate")){
        probTerminate = output["probterminate"];
      }
    }catch(...){
      cout<<"Error while parsing: Unknown output type"<<endl;
      globalIllumination = false;
      maxBounces = 1;
      probTerminate = 0.0;
      raysperpixel[0] = 4;
      raysperpixel[1] = 4;
    }
  }
}

void RayTracer::run() {
    cout<<"Running ray tracer"<<endl;    

    float width = size[0];
    float height = size[1];
    vector<double> buffer(3*width*height);
    Vector3f color(0,0,0);

    for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x++) {

        //compute direction vector for each ray
        Vector3f right = lookat.cross(up).normalized();
        float delta = 2.0 * tan((fov * M_PI / 180.0) / 2.0) / height;

        //formula given in the lectures
        Vector3f direction = lookat + tan((fov * M_PI / 180.0) / 2.0) * up - width/2.0 * delta * right + (x * delta + delta/2.0)*right - (y*delta+delta/2.0)*up;
        Ray ray(center, direction);
        color = rayCast(ray);

        buffer[3*y*width+3*x+0] = color[0];
        buffer[3*y*width+3*x+1] = color[1];
        buffer[3*y*width+3*x+2] = color[2];
      }
    }

    save_ppm(filename, buffer, width, height);
}

Vector3f RayTracer::rayCast(Ray ray){
    Vector3f intersection_point(0,0,0);
    Shape closest_shape;
    bool intersectionFound = false;
    Vector3f normal(0, 0, 0);
    Vector3f color(0, 0, 0);

    if(globalIllumination){
        // Compute the color of the ray using path tracing
        for (int s1 = 0; s1 < raysperpixel[0]; ++s1) {
          for (int s2 = 0; s2 < raysperpixel[1]; ++s2) {
              // Generate a random offset for each sample to achieve anti-aliasing
              //float random_offset_x = random_float() / (float)raysperpixel[0];
              //float random_offset_y = random_float() / (float)raysperpixel[1];
              //Ray jittered_ray = jitter_ray(ray, random_offset_x, random_offset_y);

              // Accumulate color from path tracing
              color += path_trace(ray, 0);              
          }
        }
        // Divide by the total number of samples to get the average color
        color /= (raysperpixel[0] * raysperpixel[1]);
        color = clamp(color, 0.0, 1.0);
    }
    else{
        // Compute the color of the ray using ray tracing
        intersectionFound = intersect_scene(ray, closest_shape, intersection_point, normal);
        if (intersectionFound){
          color = compute_color(closest_shape, intersection_point, normal, ray);
        }
        else{
          color = bkc;
        }
    }
    return color;
}

// Ray RayTracer::jitter_ray(const Ray& original_ray, float random_offset_x, float random_offset_y) {
//     // Assuming you have a function to generate a coordinate system based on the ray direction
//     Vector3f u, v, w;
//     coordinate_system(original_ray.d, u, v, w);

//     // Create the jittered ray
//     Vector3f jittered_direction = original_ray.d + random_offset_x * u + random_offset_y * v;
//     return Ray(original_ray.p, jittered_direction.normalized());
// }

// void RayTracer::coordinate_system(Vector3f original_ray_d, Vector3f& u, Vector3f& v, Vector3f& w) {
//     original_ray_d = original_ray_d.normalized();
//     if (std::abs(original_ray_d.x()) > std::abs(original_ray_d.y())) {
//         float inv_length = 1.0f / std::sqrt(original_ray_d.x() * original_ray_d.x() + original_ray_d.z() * original_ray_d.z());
//         u = Vector3f(-original_ray_d.z() * inv_length, 0.0f, original_ray_d.x() * inv_length);
//     } else {
//         float inv_length = 1.0f / std::sqrt(original_ray_d.y() * original_ray_d.y() + original_ray_d.z() * original_ray_d.z());
//         u = Vector3f(0.0f, original_ray_d.z() * inv_length, -original_ray_d.y() * inv_length);
//     }
//     w = original_ray_d.normalized();
//     v = u.cross(w);
// }

bool RayTracer::intersect_scene(Ray ray, Shape& hit_shape, Vector3f& intersection_pt, Vector3f& normal){
  bool intersectionFound = false;
  float t_min = numeric_limits<float>::max();
  for (const auto& sphere : spheres) {
    float t;
      if (sphereIntersection(ray, sphere, t) && t < t_min) {
        t_min = t;
        hit_shape = sphere;
        intersectionFound = true;
        normal = (ray.at(t_min) - sphere.centre);
      }
  }
  for (const auto& rectangle : rectangles) {
    float t;
      if (rectangleIntersection(ray, rectangle, t) && t < t_min) {
        t_min = t;
        hit_shape = rectangle;
        intersectionFound = true;
        normal = (rectangle.normal);
      }
  }
  if (intersectionFound){
    normal = normal.normalized();
    intersection_pt = ray.at(t_min);
  }
  return intersectionFound;
}

bool RayTracer::sphereIntersection(Ray ray, Sphere sphere, float& t) {
    Vector3f oc = ray.p - sphere.centre;
    float a = (ray.d).dot(ray.d);
    float b = 2.0 * oc.dot(ray.d);
    float c = oc.dot(oc) - sphere.radius * sphere.radius;
    float discriminant = b * b - 4 * a * c;
    if (discriminant < 0) {
        return false;
    }
    else {
        t = (-b - sqrt(discriminant)) / (2.0 * a);
        return true;
    }
}

bool RayTracer::rectangleIntersection(Ray ray, Rectangle rectangle, float& t){
  t = (rectangle.distanceToOrigin - (rectangle.normal).dot(ray.p)) / (rectangle.normal).dot(ray.d);
  Vector3f intersection_point = ray.at(t);
  if (!is_inside_triangle(intersection_point, rectangle.p1, rectangle.p2, rectangle.p3) && !is_inside_triangle(intersection_point, rectangle.p1, rectangle.p3, rectangle.p4)){
    return false;
  }
  return true;
  
}

bool RayTracer::is_inside_triangle(const Vector3f& p, const Vector3f& v0, const Vector3f& v1, const Vector3f& v2) {
    Vector3f e0 = v1 - v0;
    Vector3f e1 = v2 - v0;
    Vector3f e2 = p - v0;

    float dot00 = e0.dot(e0);
    float dot01 = e0.dot(e1);
    float dot02 = e0.dot(e2);
    float dot11 = e1.dot(e1);
    float dot12 = e1.dot(e2);

    float inv_denom = 1.0 / (dot00 * dot11 - dot01 * dot01);
    float u = (dot11 * dot02 - dot01 * dot12) * inv_denom;
    float v = (dot00 * dot12 - dot01 * dot02) * inv_denom;

    float epsilon = 1e-4;
    return (u >= -epsilon && v >= -epsilon && u + v <= 1.0 + epsilon);
}

Vector3f RayTracer::compute_color(Shape shape, Vector3f intersection_pt, Vector3f normal, Ray ray) {
    Vector3f color(0.0, 0.0, 0.0);
    Vector3f ambient(0.0, 0.0, 0.0);
    Vector3f diffuse(0.0, 0.0, 0.0);
    Vector3f specular(0.0, 0.0, 0.0);
    // Compute the ambient color
    ambient = shape.ka * shape.ac;

    for (PointLight light : pointLights) {
      Vector3f L = light.centre - intersection_pt; 
      float max_distance = L.norm();
      L = L.normalized();
      Ray ptL(intersection_pt + L * 1e-3, L);
      if (normal.dot(L) > 0 && !testShadow(ptL, max_distance)) {
          // Calculate diffuse illumination
          diffuse += (shape.dc).cwiseProduct(light.id * normal.dot(L))* shape.kd;
          // Calculate specular illumination
          Vector3f Cv = ray.p - intersection_pt; 
          Cv = Cv.normalized();
          //Vector3f Rv = (2 * normal * normal.dot(L)) - L;
          //Rv = Rv.normalized();
          Vector3f H = (L+Cv)/(L+Cv).norm();
          if (H.dot(normal) > 0) {
              specular += (shape.sc).cwiseProduct(light.is * pow(H.dot(normal), shape.pc))* shape.ks;
          }
          // if (Cv.dot(Rv) > 0) {
          //   specular += (shape.sc).cwiseProduct(light.is * pow(Cv.dot(Rv), shape.pc))* shape.ks;
          // }
                  
      }
    }

    for (AreaLight light : areaLights) {
      if(light.useCenter || globalIllumination){
        light.GetCenter();
        Vector3f L = light.centre - intersection_pt; 
        float max_distance = L.norm();
        L = L.normalized();
        Ray ptL(intersection_pt + L * 1e-3, L);
        if (normal.dot(L) > 0 && !testShadow(ptL, max_distance)) {
          // Calculate diffuse illumination
          diffuse += (shape.dc).cwiseProduct(light.id * normal.dot(L))* shape.kd;
          // Calculate specular illumination
          Vector3f Cv = ray.p - intersection_pt; 
          Cv = Cv.normalized();
          Vector3f Rv = (2 * normal * normal.dot(L)) - L;
          Rv = Rv.normalized();
          Vector3f H = (L+Cv)/(L+Cv).norm();
          if (Cv.dot(Rv) > 0) {
              specular += (shape.sc).cwiseProduct(light.is * pow(Cv.dot(Rv), shape.pc))* shape.ks;
          }
          // if (H.dot(normal) > 0) {
          //     specular += (shape.sc).cwiseProduct(light.is * pow(H.dot(normal), shape.pc))* shape.ks;
          // }          
        }
      }
      else{
        int a = raysperpixel[0];
        int b = raysperpixel[1];

        for (int i = 0; i < a; ++i) {
          for (int j = 0; j < a; ++j) {
            for (int k = 0; k < b; ++k) {
              float u = (i + random_float()) / (float)a;
              float v = (j + random_float()) / (float)a;
              
              Vector3f p_uv = light.p1 * (1 - u) * (1 - v) + light.p2 * u * (1 - v) + light.p3 * u * v + light.p4 * (1 - u) * v;
              Vector3f L = p_uv - intersection_pt;
              float max_distance = L.norm();
              L = L.normalized();
              
              Ray ptL(intersection_pt + L * 1e-3, L);
              if (normal.dot(L) > 0 && !testShadow(ptL, max_distance)) {
                diffuse += (shape.dc).cwiseProduct(light.id * normal.dot(L)) * shape.kd / (a * a * b);

                Vector3f Cv = ray.p - intersection_pt;
                Cv = Cv.normalized();
                Vector3f H = (L + Cv) / (L + Cv).norm();
                if (H.dot(normal) > 0) {
                  specular += (shape.sc).cwiseProduct(light.is * pow(H.dot(normal), shape.pc)) * shape.ks / (a * a * b);
                }
              }
            }
          }
        }
      }
    }
    color = ambient + diffuse + specular;
    
    // Clamp the color values to [0, 1]
    color = clamp(color, 0.0, 1.0);
    
    return color;
}

Vector3f RayTracer::path_trace(Ray& ray, int depth) {
    if (depth > maxBounces) {
        return Vector3f(0.0, 0.0, 0.0);
    }

    Shape hit_shape;
    Vector3f intersection_pt(0,0,0);
    Vector3f normal(0,0,0);
    if (!intersect_scene(ray, hit_shape, intersection_pt, normal)) {
        return bkc;
    }

    //Russian Roulette termination
    if (random_float() < probTerminate) {
        return hit_shape.ac * hit_shape.ka;
    }

    // Calculate direct illumination from point/area lights (diffuse only)
    Vector3f direct_illumination(0.0, 0.0, 0.0);
    direct_illumination = compute_color(hit_shape, intersection_pt, normal, ray);

    //Sample a random direction for the next ray bounce
    Vector3f wi = sample_hemisphere(normal);

    //Compute the BRDF
    Vector3f brdf = hit_shape.kd * hit_shape.dc / M_PI;

    //Create the next ray
    Ray next_ray(intersection_pt + wi * 1e-3, wi);

    //Compute the cosine term
    float cos_theta = std::max(normal.dot(wi), 0.0f);

    //Calculate the path contribution using Monte Carlo integration
    Vector3f path_contribution = path_trace(next_ray, depth + 1);

    return hit_shape.ac * hit_shape.ka + direct_illumination + (1.0f - probTerminate) * brdf.cwiseProduct(path_contribution) * cos_theta;
    //return direct_illumination;
}

Vector3f RayTracer::clamp(Vector3f v, float min, float max) {
    for (int i = 0; i < 3; i++) {
        if (v[i] < min) {
            v[i] = min;
        }
        else if (v[i] > max) {
            v[i] = max;
        }
    }
    return v;
}

// Returns a random float between 0 and 1st
float RayTracer::random_float() {
    static std::random_device rd;
    static std::mt19937 generator(rd());
    static std::uniform_real_distribution<float> distribution(0.0, 1.0);
    return distribution(generator);
}

bool RayTracer::testShadow(Ray &ray, float max_distance) {

    // Iterate over every face in every sphere and rectangle to test for ray intersection
    float t = 0.0;
    for (const auto& sphere : spheres) {
        if (sphereIntersection(ray, sphere, t)) {
          if(t>0.00001 && t<max_distance)
            return true;
        }
    }
    for (const auto& rectangle : rectangles) {
        if (rectangleIntersection(ray, rectangle, t)) {
          if(t>0.00001 && t<max_distance)
            return true;
        }
    }

    // Return false if no intersection was found
    return false;
}

Vector3f RayTracer::sample_hemisphere(const Vector3f& normal) {
    float u = random_float();
    float v = random_float();

    float theta = 2 * M_PI * u;
    float phi = acos(2 * v - 1);

    Vector3f direction(sinf(phi) * cosf(theta), sinf(phi) * sinf(theta), cosf(phi));

    // Create an orthonormal basis with the normal as one of the axes
    Vector3f w = normal.normalized();
    Vector3f u_temp = (fabs(w.x()) > 0.1) ? Vector3f(0, 1, 0) : Vector3f(1, 0, 0);
    Vector3f v_temp = w.cross(u_temp);
    Vector3f u_axis = w.cross(v_temp).normalized();
    Vector3f v_axis = w.cross(u_axis);

    // Convert the direction to the orthonormal basis
    Vector3f sampled_direction = direction.x() * u_axis + direction.y() * v_axis + direction.z() * w;

    return sampled_direction.normalized();
}