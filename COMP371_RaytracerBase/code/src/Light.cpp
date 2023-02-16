#include "Light.h"

void Light::SetLightIntensity(float id1, float id2, float id3, float is1, float is2, float is3) {
    this->id[0] = id1;
    this->id[1] = id2;
    this->id[2] = id3;
    this->is[0] = is1;
    this->is[1] = is2;
    this->is[2] = is3;
}