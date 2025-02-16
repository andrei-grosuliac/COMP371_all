#include "Sphere.h"
#include "Rectangle.h"
#include "RayTracer.h"
#include "AreaLight.h"
#include "PointLight.h"
#include <random>
//#include <thread>
//#include <pthread.h>

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
      if(light.contains("use")) {
        if(light["use"] == false) {
          continue;
        }
      }
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
  for (auto& file : scene["output"]) {
      File f;
      f.filename = file["filename"];
      f.size[0] = file["size"][0];
      f.size[1] = file["size"][1];

      f.lookat = Vector3f(file["lookat"][0], file["lookat"][1], file["lookat"][2]);
      f.up = Vector3f(file["up"][0], file["up"][1], file["up"][2]);
      f.center = Vector3f(file["centre"][0], file["centre"][1], file["centre"][2]);
      f.bkc = Vector3f(file["bkc"][0], file["bkc"][1], file["bkc"][2]);
      f.fov = file["fov"];
      try{
        if(file.contains("raysperpixel")){
          f.raysperpixel[0] = file["raysperpixel"][0];
          f.raysperpixel[1] = file["raysperpixel"][1];
        }
        else{
          f.raysperpixel[0] = 5;
          f.raysperpixel[1] = 5;
        }
        if(file.contains("globalillum")){
          f.globalIllumination = file["globalillum"];
        }
        else{
          f.globalIllumination = false;
        }
        if(file.contains("maxbounces")){
          f.maxBounces = file["maxbounces"];
        }
        if(file.contains("probterminate")){
          f.probTerminate = file["probterminate"];
        }
        if(file.contains("antialiasing")){
          f.antiAliasing = file["antialiasing"];
        }
        else{
          f.antiAliasing = false;
        }
      }catch(...){
        cout<<"Error while parsing: Unknown output type"<<endl;
        f.globalIllumination = false;
        f.maxBounces = 1;
        f.antiAliasing = false;
        f.probTerminate = 0.0;
        f.raysperpixel[0] = 5;
        f.raysperpixel[1] = 5;        
      }
      files.push_back(f);        
  }
}

void RayTracer::run() {
    for (auto& file : files) {
      checkAreaLights(file);
        cout<<"Running ray tracer"<<endl;    

        float width = file.size[0];
        float height = file.size[1];
        vector<double> buffer(3*width*height);

        //Uncomment to use multithreading

        //int numThreads = std::thread::hardware_concurrency(); // Get the number of supported threads
        //numThreads = 1; // Get the number of supported threads
        //std::vector<std::thread> threads(numThreads);         // Create a vector of threads

        //int sectionHeight = height / numThreads; // Divide the image into sections

        // Create and launch threads
        //for (int i = 0; i < numThreads; i++) {
            //int startY = i * sectionHeight;
           // int endY = (i == numThreads - 1) ? height : (i + 1) * sectionHeight;
           // threads[i] = std::thread(&RayTracer::render_section, this, startY, endY, width, height, std::ref(buffer),  std::ref(file));
        //}

        // Join the threads (wait for them to finish)
        //for (auto& thread : threads) {
        //    thread.join();
        //}

        render_section(0, 0, width, height, buffer, file);

        save_ppm(file.filename, buffer, width, height);
    }
    
}

void RayTracer::render_section(int startY, int endY, int imageWidth, int imageHeight, vector<double>& buffer, File& file) {
    int numSamples = file.raysperpixel[0];
    int raysPerCell = file.raysperpixel[1];
    for (int y = 0; y < imageHeight; y++) {
        for (int x = 0; x < imageWidth; x++) {          
          Vector3f accumulated_color(0, 0, 0);

          if(file.antiAliasing || file.globalIllumination){
            for (int sy = 0; sy < numSamples; sy++) {
              for (int sx = 0; sx < numSamples; sx++) {
                for (int i = 0; i < raysPerCell; ++i) {
                  // Stratify the pixel by dividing it into a grid
                  float u = (sx + random_float()) / numSamples;
                  float v = (sy + random_float()) / numSamples;

                  // Compute direction vector for each ray
                  Vector3f right = file.lookat.cross(file.up).normalized();
                  float delta = 2.0 * tan((file.fov * M_PI / 180.0) / 2.0) / imageHeight;

                  // Formula given in the lectures, modified to use stratified sampling
                  Vector3f direction = file.lookat + tan((file.fov * M_PI / 180.0) / 2.0) * file.up
                      - imageWidth / 2.0 * delta * right
                      + ((x + u) * delta) * right
                      - ((y + v) * delta) * file.up;
                  direction = direction.normalized();

                  Ray ray(file.center, direction);
                  Vector3f color = rayCast(ray, file);

                  accumulated_color += color;
                }
              }
            }
            // Average the accumulated color
            accumulated_color /= (numSamples * numSamples * raysPerCell);
          }
          else{
            // Compute direction vector for each ray
            Vector3f right = file.lookat.cross(file.up).normalized();
            float delta = 2.0 * tan((file.fov * M_PI / 180.0) / 2.0) / imageHeight;

            // Formula given in the lectures, modified to use stratified sampling
            Vector3f direction = file.lookat + tan((file.fov * M_PI / 180.0) / 2.0) * file.up
                - imageWidth / 2.0 * delta * right
                + (x * delta) * right
                - (y * delta) * file.up;
            direction = direction.normalized();

            Ray ray(file.center, direction);
            accumulated_color = rayCast(ray, file);
          }            

          buffer[3 * y * imageWidth + 3 * x + 0] = accumulated_color[0];
          buffer[3 * y * imageWidth + 3 * x + 1] = accumulated_color[1];
          buffer[3 * y * imageWidth + 3 * x + 2] = accumulated_color[2];
        }
    }
}

void RayTracer::checkAreaLights(File& file){
  for (auto& light: areaLights){
    if (!light.useCenter){
      file.antiAliasing = false;
    }
  }
}

Vector3f RayTracer::rayCast(Ray& ray, File& file){
    Vector3f intersection_point(0,0,0);
    Shape closest_shape;
    bool intersectionFound = false;
    Vector3f normal(0, 0, 0);
    Vector3f color(0, 0, 0);

    if(file.globalIllumination){
        color = path_trace(ray, 0, file);
        color = clamp(color, 0.0, 1.0);        
    }
    else{
        // Compute the color of the ray using ray tracing
        intersectionFound = intersect_scene(ray, closest_shape, intersection_point, normal);
        if (intersectionFound){
          color = compute_color(closest_shape, intersection_point, normal, ray, file);
        }
        else{
          color = file.bkc;
        }
    }
    return color;
}

bool RayTracer::intersect_scene(Ray& ray, Shape& hit_shape, Vector3f& intersection_pt, Vector3f& normal){
  bool intersectionFound = false;
  float t_min = numeric_limits<float>::max();
  for (auto& sphere : spheres) {
    float t;
      if (sphereIntersection(ray, sphere, t) && t < t_min) {
        t_min = t;
        hit_shape = sphere;
        intersectionFound = true;
        normal = (ray.at(t_min) - sphere.centre).normalized();
      }
  }
  for (auto& rectangle : rectangles) {
    float t;
      if (rectangleIntersection(ray, rectangle, t) && t < t_min) {
        t_min = t;
        hit_shape = rectangle;
        intersectionFound = true;
        normal = (rectangle.normal);
      }
  }
  if (intersectionFound){
    intersection_pt = ray.at(t_min);
  }
  return intersectionFound;
}

bool RayTracer::sphereIntersection(Ray& ray, Sphere& sphere, float& t) {
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
		if(t<0.001){
		  return false;
		}
        return true;
    }
}

bool RayTracer::rectangleIntersection(Ray& ray, Rectangle& rectangle, float& t) {
    // Check if ray is parallel to the rectangle plane
    float ndotd = rectangle.normal.dot(ray.d);
    if (std::abs(ndotd) < 1e-6) {
        return false;
    }

    // Calculate the intersection point on the plane
    t = (rectangle.normal.dot(rectangle.p1) - rectangle.normal.dot(ray.p)) / ndotd;
    if(t<0.001){
      return false;
    }
    Vector3f intersection_point = ray.at(t);


    // Check if intersection point is inside the rectangle
    Vector3f edge1 = rectangle.p2 - rectangle.p1;
    Vector3f edge2 = rectangle.p4 - rectangle.p1;
    Vector3f V = intersection_point - rectangle.p1;

    float d00 = edge1.dot(edge1);
    float d01 = edge1.dot(edge2);
    float d02 = edge1.dot(V);
    float d11 = edge2.dot(edge2);
    float d12 = edge2.dot(V);

    float invDenom = 1.0 / (d00 * d11 - d01 * d01);
    float u = (d11 * d02 - d01 * d12) * invDenom;
    float v = (d00 * d12 - d01 * d02) * invDenom;

    if (u >= 0.0f && u <= 1.0f && v >= 0.0f && v <= 1.0f) {
        return true;
    }

    return false;
}

Vector3f RayTracer::compute_color(Shape& shape, Vector3f& intersection_pt, Vector3f& normal, Ray& ray, File& file) {
    Vector3f color(0.0, 0.0, 0.0);    
    Vector3f diffuse(0.0, 0.0, 0.0);
    Vector3f specular(0.0, 0.0, 0.0);
    Vector3f ambient(0.0, 0.0, 0.0);

    // Compute the ambient color
    if(!file.globalIllumination){
      ambient = shape.ka * shape.ac;
    }

    for (PointLight light : pointLights) {
      Vector3f L = light.centre - intersection_pt; 
      float max_distance = L.norm();
      L = L.normalized();
      Ray ptL(intersection_pt + L * 1e-3, L);
      if (normal.dot(L) > 0 && !testShadow(ptL, max_distance)) {
          // Calculate diffuse illumination
          diffuse += (shape.dc).cwiseProduct(light.id * normal.dot(L))* shape.kd;

          // Calculate specular illumination
          if(!file.globalIllumination){            
            Vector3f Cv = ray.p - intersection_pt; 
            Cv = Cv.normalized();
            Vector3f H = (L+Cv)/(L+Cv).norm();
            if (H.dot(normal) > 0) {
                specular += (shape.sc).cwiseProduct(light.is * pow(H.dot(normal), shape.pc))* shape.ks;
            }
          }
                  
      }
    }

    for (AreaLight light : areaLights) {
      if(light.useCenter || file.globalIllumination){

        light.GetCenter();
        Vector3f L = light.centre - intersection_pt; 
        float max_distance = L.norm();
        L = L.normalized();
        Ray ptL(intersection_pt + L * 1e-3, L);

        if (normal.dot(L) > 0 && !testShadow(ptL, max_distance)) {

          // Calculate diffuse illumination
          diffuse += (shape.dc).cwiseProduct(light.id * normal.dot(L))* shape.kd;

          if(!file.globalIllumination){
            // Calculate specular illumination
            Vector3f Cv = ray.p - intersection_pt; 
            Cv = Cv.normalized();
            Vector3f Rv = (2 * normal * normal.dot(L)) - L;
            Rv = Rv.normalized();
            Vector3f H = (L+Cv)/(L+Cv).norm();
            if (Cv.dot(Rv) > 0) {
                specular += (shape.sc).cwiseProduct(light.is * pow(Cv.dot(Rv), shape.pc))* shape.ks;
            } 
          }      
        }
      }
      else{
        int a = file.raysperpixel[0];
        int b = file.raysperpixel[1];

        for (int i = 0; i < a; ++i) {
          for (int i = 0; i < a; ++i) {
            for (int j = 0; j < b; ++j) {
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
    if(!file.globalIllumination){
    color = ambient + diffuse + specular;
    }
    else{
      color = diffuse;
    }
    
    // Clamp the color values to [0, 1]
    color = clamp(color, 0.0, 1.0);
    
    return color;
}

Vector3f RayTracer::path_trace(Ray& ray, int depth, File& file) {
    if (depth >= file.maxBounces) {
        return Vector3f(0.0, 0.0, 0.0);
    }
   
    Shape hit_shape;
    Vector3f intersection_pt(0,0,0);
    Vector3f normal(0,0,0);
    if (!intersect_scene(ray, hit_shape, intersection_pt, normal)) {
      return file.bkc;
    }
    
    Vector3f direct_illumination(0.0, 0.0, 0.0);
    Vector3f indirect_illumination(0.0, 0.0, 0.0);

    //Calculate direct illumination from point/area lights (diffuse only)
    direct_illumination = compute_color(hit_shape, intersection_pt, normal, ray, file); 
    
    //Russian Roulette termination
    if (random_float() < file.probTerminate) {      
      return direct_illumination;
    }
    //Sample a random direction for the next ray bounce
    Vector3f wi = sample_hemisphere(normal);
    
    //Compute the cosine term
    float cos_theta = std::max(normal.dot(wi), 0.0f);

    //Compute brdf
    Vector3f brdf = hit_shape.dc * hit_shape.kd / M_PI;

    //Create the next ray
    Ray next_ray(intersection_pt + wi * 1e-3, wi);
    indirect_illumination += ((brdf.cwiseProduct(path_trace(next_ray, depth + 1, file)) * cos_theta * (2*M_PI)));
    
    return (direct_illumination + indirect_illumination);
}

Vector3f RayTracer::clamp(Vector3f& v, float min, float max) {
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
    for (auto& sphere : spheres) {
        if (sphereIntersection(ray, sphere, t)) {
          if(t>0.00001 && t<max_distance)
            return true;
        }
    }
    for (auto& rectangle : rectangles) {
        if (rectangleIntersection(ray, rectangle, t)) {
          if(t>0.00001 && t<max_distance)
            return true;
        }
    }

    // Return false if no intersection was found
    return false;
}

Vector3f RayTracer::sample_hemisphere(Vector3f& normal) {
    // Generate two random numbers between 0 and 1
    float u1 = random_float();
    float u2 = random_float();

    // Generate a random direction on the unit disk (x, y plane)
    float r = sqrt(u1);
    float theta = 2 * M_PI * u2;
    float x = r * cos(theta);
    float y = r * sin(theta);

    // Compute the z component
    float z = sqrt(1 - x * x - y * y);

    // Create a local coordinate system with the normal as the z-axis
    Vector3f tangent = (fabs(normal.x()) > 0.1) ? Vector3f(0, 1, 0) : Vector3f(1, 0, 0);
    Vector3f bitangent = normal.cross(tangent).normalized();
    tangent = bitangent.cross(normal);

    // Convert the direction from the local coordinate system to the global coordinate system
    Vector3f sampled_direction = x * tangent + y * bitangent + z * normal;

    return sampled_direction.normalized();
}


int RayTracer::testHemisphere() {
    Vector3f normal(0, 1, 0); // example normal
    int num_samples = 1000;
    int valid_samples = 0;

    for (int i = 0; i < num_samples; ++i) {
        Vector3f sample = sample_hemisphere(normal);
        if (sample.dot(normal) >= 0.0) {
            valid_samples++;
        } else {
            std::cout << "Invalid sample: " << sample << std::endl;
        }
    }

    std::cout << "Valid samples: " << valid_samples << " / " << num_samples << std::endl;

    return 0;
}