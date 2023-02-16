#include "Sphere.h"
#include "Rectangle.h"
#include "RayTracer.h"
#include "AreaLight.h"
#include "PointLight.h"

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

            lights.push_back(pointLight);
        }
        else if (light["type"] == "area") {
            cout<<"Creating area light"<<endl;

            AreaLight areaLight(light["p1"][0], light["p1"][1], light["p1"][2], 
            light["p2"][0], light["p2"][1], light["p2"][2], 
            light["p3"][0], light["p3"][1], light["p3"][2], 
            light["p4"][0], light["p4"][1], light["p4"][2]);

            areaLight.SetLightIntensity(light["id"][0], light["id"][1], light["id"][2], light["is"][0], light["is"][1], light["is"][2]);

            lights.push_back(areaLight);
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
    fov = output["fov"];
  }
}

void RayTracer::run() {
    cout<<"Running ray tracer"<<endl;    

    float width = size[0];
    float height = size[1];
    vector<double> buffer(3*width*height);

    for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x++) {

        //compute direction vector for each ray
        Vector3f right = lookat.cross(up).normalized();
        float delta = 2.0 * tan((fov * M_PI / 180.0) / 2.0) / height;
        //formula given in the lectures
        Vector3f direction = center + lookat + tan((fov * M_PI / 180.0) / 2.0) * up - width/2.0 * delta * right + (x * delta + delta/2.0)*right - (y*delta+delta/2.0)*up;
        //cout<<direction[0]<<", "<< direction[1]<<", "<< direction[2]<<endl;

        Ray ray(center, direction);
        rayCast(ray);

        buffer[3*y*width+3*x+0] = color[0];
        buffer[3*y*width+3*x+1] = color[1];
        buffer[3*y*width+3*x+2] = color[2];
      }
    }

    save_ppm(filename, buffer, width, height);
}

void RayTracer::rayCast(Ray ray){
    float t_min = numeric_limits<float>::max();
    Vector3f closest_color = bkc;
    for (const auto& sphere : spheres) {
      float t;
        if (sphereIntersection(ray, sphere, t) && t < t_min) {
        t_min = t;
        closest_color = {1,1,1};
        }
    }  
    color = closest_color;
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