#include "PointLight.h"

// Constructor
PointLight::PointLight(float x, float y, float z) {
    cout << "Created point light at (" << x << ", " << y << ", " << z << ")." << endl;
    centre[0] = x;
    centre[1] = y;
    centre[2] = z;
}