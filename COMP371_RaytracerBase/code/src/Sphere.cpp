#include "Sphere.h"


// Constructor
Sphere::Sphere(float x, float y, float z, float radius) {
    cout << "Created sphere with radius " << radius << " at (" << x << ", " << y << ", " << z << ")." << endl;
    centre[0] = x;
    centre[1] = y;
    centre[2] = z;
    this->radius = radius;
    this->type = "sphere";
}