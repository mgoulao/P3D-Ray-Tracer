#ifndef CAMERA_H
#define CAMERA_H

#include <cmath>
#include <stdio.h>
using namespace std;

#include "vector.h"
#include "ray.h"


#define PI				3.141592653589793238462f

class Camera
{

private:
	Vector eye, at, up; 
	float fovy, vnear, vfar, plane_dist, focal_ratio, aperture;
	float w, h;
	int res_x, res_y;
	Vector u, v, n;

public:
	Vector GetEye() { return eye; }
	int GetResX()  { return res_x; }
    int GetResY()  { return res_y; }
	float GetFov() { return fovy; }
	float GetPlaneDist() { return plane_dist; }
	float GetFar() {return vfar; }
	float GetAperture() { return aperture; }

    Camera( Vector from, Vector At, Vector Up, float angle, float hither, float yon, int ResX, int ResY, float Aperture_ratio, float Focal_ratio) {
	    eye = from;
	    at = At;
	    up = Up;
	    fovy = angle;
	    vnear = hither;
	    vfar = yon;
	    res_x = ResX;
	    res_y = ResY;
		focal_ratio = Focal_ratio;

        // set the camera frame uvn
        n = ( eye - at );
        plane_dist = n.length();
	    n = n / plane_dist;

	    u = up % n;
	    u = u / u.length();

	    v = n % u;

        //Dimensions of the vis window
	    h = 2 * plane_dist * tan( (PI * angle / 180) / 2.0f );
        w = ( (float) res_x / res_y ) * h;  

		aperture = Aperture_ratio * (w / res_x); //Lens aperture = aperture_ratio * pixel_size

		printf("\nwidth=%f height=%f fov=%f, viewplane distance=%f, pixel size=%.3f\n", w,h, fovy,plane_dist, w/res_x);
		if (Aperture_ratio != 0) printf("\nDepth-Of-Field effect enabled with a lens aperture = %.1f\n", Aperture_ratio);
    }

	void SetEye(Vector from) {
		eye = from;
		// set the camera frame uvn
		n = (eye - at);
		plane_dist = n.length();
		n = n / plane_dist;
		u = up % n;
		u = u / u.length();
		v = n % u;
	}

	Ray PrimaryRay(const Vector& pixel_sample) //  Rays cast from the Eye to a pixel sample which is in Viewport coordinates
	{
		Vector ray_dir;

		float r =  w / 2;
		float t =  h / 2;
		float l = -w / 2;
		float b = -h / 2;

		float u_s = l + (r - l) * (pixel_sample.x/this->GetResX());
		float v_s = b + (t - b) * (pixel_sample.y/this->GetResY());
		float n_s = -(eye - at).length();

		Vector s = transforUVWToXYZ(Vector(u_s, v_s, n_s));

		ray_dir = s - eye;
		ray_dir.normalize();
		
		return Ray(eye, ray_dir);  
	}

	Ray PrimaryRay(const Vector& lens_sample, const Vector& pixel_sample) // DOF: Rays cast from  a thin lens sample to a pixel sample
	{
		Vector rayDir;
		Vector eyeOffset;
		Vector pixel = getPixelCoordinates(pixel_sample);

		float focalDistance = plane_dist * focal_ratio;
		Vector p = Vector(pixel.x * (focalDistance / plane_dist), pixel.y * (focalDistance / plane_dist), -focalDistance);

		Vector lens = lens_sample;
		eyeOffset = transforUVWToXYZ(lens*aperture);
		rayDir = (transforUVWToXYZ(p) - eyeOffset).normalize();
		return Ray(eyeOffset, rayDir);
	}

	Vector getPixelCoordinates(const Vector& pixelSample) {
		float r = w / 2;
		float t = h / 2;
		float l = -w / 2;
		float b = -h / 2;

		float u_s = l + (r - l) * (pixelSample.x / GetResX());
		float v_s = b + (t - b) * (pixelSample.y / GetResY());
		float n_s = -(eye - at).length();

		return Vector(u_s, v_s, n_s);
	}

	Vector transforUVWToXYZ(Vector uvw) {
		Vector s; // transform s from uvw to xyz
		s.x = eye.x + uvw.x * u.x + uvw.y * v.x + uvw.z * n.x;
		s.y = eye.y + uvw.x * u.y + uvw.y * v.y + uvw.z * n.y;
		s.z = eye.z + uvw.x * u.z + uvw.y * v.z + uvw.z * n.z;
		return s;
	}

};

#endif