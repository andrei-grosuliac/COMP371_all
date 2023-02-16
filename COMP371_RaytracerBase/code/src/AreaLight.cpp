#include "AreaLight.h"

// Constructor
AreaLight::AreaLight(float p1a, float p1b, float p1c, float p2a, float p2b, float p2c, float p3a, float p3b, float p3c, float p4a, float p4b, float p4c) {
    cout << "Created area light with vertices (" << p1a << ", " << p1b << ", " << p1c << ")," 
    "(" << p2a << ", " << p2b << ", " << p2c << ")," 
    "(" << p3a << ", " << p3b << ", " << p3c << ")," 
    "(" << p4a << ", " << p4b << ", " << p4c << ")." << endl;
    
    p1[0] = p1a;
    p1[1] = p1b;
    p1[2] = p1c;
    p2[0] = p2a;
    p2[1] = p2b;
    p2[2] = p2c;
    p3[0] = p3a;
    p3[1] = p3b;
    p3[2] = p3c;
    p4[0] = p4a;
    p4[1] = p4b;
    p4[2] = p4c;
}