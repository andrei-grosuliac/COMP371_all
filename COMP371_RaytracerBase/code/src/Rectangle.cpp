#include "Rectangle.h"


// Constructor
Rectangle::Rectangle(float p1a, float p1b, float p1c, float p2a, float p2b, float p2c, float p3a, float p3b, float p3c, float p4a, float p4b, float p4c) {
    cout << "Created rectangle with points" 
    "(" << p1a << ", " << p1b << ", " << p1c << ")," 
    "(" << p2a << ", " << p2b << ", " << p2c << ")," 
    "(" << p3a << ", " << p3b << ", " << p3c << ")," 
    "(" << p4a << ", " << p4b << ", " << p4c << ")." << endl;

    this->p1[0] = p1a;
    this->p1[1] = p1b;
    this->p1[2] = p1c;
    this->p2[0] = p2a;
    this->p2[1] = p2b;
    this->p2[2] = p2c;
    this->p3[0] = p3a;
    this->p3[1] = p3b;
    this->p3[2] = p3c;
    this->p4[0] = p4a;
    this->p4[1] = p4b;
    this->p4[2] = p4c;

    this->type = "rectangle";
}