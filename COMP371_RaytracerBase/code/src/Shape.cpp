#include "Shape.h"

void Shape::SetReflectionCoefficients(float ka, float kd, float ks) {
    this->ka = ka;
    this->kd = kd;
    this->ks = ks;
}

void Shape::SetReflectionColor(float ac1, float ac2, float ac3, float dc1, float dc2, float dc3, float sc1, float sc2, float sc3) {
    this->ac[0] = ac1;
    this->ac[1] = ac2;
    this->ac[2] = ac3;
    this->dc[0] = dc1;
    this->dc[1] = dc2;
    this->dc[2] = dc3;
    this->sc[0] = sc1;
    this->sc[1] = sc2;
    this->sc[2] = sc3;
}

void Shape::SetPhongCoefficients(float pc) {
    this->pc = pc;
}