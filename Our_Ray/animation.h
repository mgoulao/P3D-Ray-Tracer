#ifndef ANIMATION_H
#define ANIMATION_H

class Sphere; // forward declaration

class Animation {
private:
	Vector initialPosition;
	Vector finalPosition;
	Vector stepVec;
	int frames;
	Sphere* object;
public:
	Animation(Vector& initialPosition, Vector& finalPosition, int frames);

	void step(int frame);
	void addObject(Sphere* obj);
};

#endif