#include "vector.h"
#include "scene.h"
#include "animation.h"

Animation::Animation(Vector& initialPosition, Vector& finalPosition, int frames) {
	this->finalPosition = finalPosition;
	this->initialPosition = initialPosition;
	this->frames = frames;
	this->stepVec = (finalPosition - initialPosition) / frames;
}

void Animation::step(int frame) {
	if (object != NULL) {
		Vector currentPosition = initialPosition + stepVec * frame;
		object->translate(currentPosition);
	}
}

void Animation::addObject(Sphere* obj) {
	object = obj;
	step(0); // override object position
}