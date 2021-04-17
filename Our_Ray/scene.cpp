#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <IL/il.h>
#include <math.h>
#include "maths.h"
#include "scene.h"
#include "sampler.h"
#include "rayAccelerator.h"

Vector Light::sampleLight(Vector passSample) {
	return position;
}

Vector AreaLight::sampleLight(Vector passSample) {
	P0;
	P1;
	return position + P0 * (passSample.x - 0.5) + P1 * (passSample.y - 0.5);
}

vector<Light*> AreaLight::decompose(int n) {
	float factor = 1 / (float) n;
	vector<Light*> decomposedLights;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			Vector point = P0 * (i + 0.5) * factor + P1 * (j + 0.5) * factor + position;
			decomposedLights.push_back(new Light(point, color/(n*n)));
		}
	}
	return decomposedLights;
}

Triangle::Triangle(Vector& P0, Vector& P1, Vector& P2)
{
	points[0] = P0; points[1] = P1; points[2] = P2;

	/* Calculate the normal */
	Vector P01 = P1 - P0;
	Vector P02 = P2 - P0;
	normal = P01 % P02;
	normal.normalize();

	//Calculate the Min and Max for bounding box
	float minX = min(P0.x, min(P1.x, P2.x));
	float minY = min(P0.y, min(P1.y, P2.y));
	float minZ = min(P0.z, min(P1.z, P2.z));
	float maxX = max(P0.x, max(P1.x, P2.x));
	float maxY = max(P0.y, max(P1.y, P2.y));
	float maxZ = max(P0.z, max(P1.z, P2.z));
	Min = Vector(minX, minY, minZ);
	Max = Vector(maxX, maxY, maxZ);


	// enlarge the bounding box a bit just in case...
	Min -= EPSILON;
	Max += EPSILON;
}

AABB Triangle::GetBoundingBox() {
	return(AABB(Min, Max));
}

Vector Triangle::getNormal(Vector point)
{
	return normal;
}

//
// Ray/Triangle intersection test using Tomas Moller-Ben Trumbore algorithm.
//

bool Triangle::intercepts(Ray& r, float& t) {
	float a = points[0].x - points[1].x;
	float b = points[0].y - points[1].y;
	float c = points[0].z - points[1].z;
	float d = points[0].x - points[2].x;
	float e = points[0].y - points[2].y;
	float f = points[0].z - points[2].z;
	float g = r.direction.x;
	float h = r.direction.y;
	float i = r.direction.z;
	float j = points[0].x - r.origin.x;
	float k = points[0].y - r.origin.y;
	float l = points[0].z - r.origin.z;
	float M = a * (e * i - h * f) + b* (g * f - d * i) + c * (d * h - e * g);
	t = - (f * (a * k - j * b) + e * (j * c - a * l) + d * (b * l - k * c)) / M;
	if (t < 0) {
		return false;
	}
	float gamma = (i * (a * k - j * b) + h * (j * c - a * l) + g * (b * l - k * c)) / M;
	if (gamma < 0 || gamma > 1) {
		return false;
	}
	float beta = (j * (e * i - h * f) + k * (g * f - d * i) + l * (d * h - e * g)) / M;
	if (beta < 0 || beta > 1 - gamma) {
		return false;
	}
	return true;
}

Plane::Plane(Vector& a_PN, float a_D)
	: PN(a_PN), D(a_D)
{}

Plane::Plane(Vector& P0, Vector& P1, Vector& P2)
{
   float l;

   //Calculate the normal plane: counter-clockwise vectorial product.
   Vector P01 = P1 - P0;
   Vector P02 = P2 - P0;
   PN = P01 % P02;		

   if ((l=PN.length()) == 0.0)
   {
     cerr << "DEGENERATED PLANE!\n";
   }
   else
   {
     PN.normalize();
	 P = P0;
	 //Calculate D
     D = 0; // TODO
   }
}

//
// Ray/Plane intersection test.
//

bool Plane::intercepts( Ray& r, float& t )
{
	float denominator = r.direction * PN;
	
	if (denominator > 0) return false;
	t = ((P - r.origin) * PN) / denominator;
	
	return t > 0;
}

Vector Plane::getNormal(Vector point) 
{
	return PN;
}


bool Sphere::intercepts(Ray& r, float& t)
{
	Vector d = r.direction;
	Vector e = r.origin;
	Vector c = this->center;
	float R = this->radius;
	t = 0; 

	//if ((c - e) * d < 0) { // determine first if the sphere is behind the ray origin
	//	return false;
	//}

	float discriminator = pow((d * (e - c)), 2) - (d * d) * ((e - c) * (e - c) - R * R);
	if (discriminator < 0) return false;
	float t_minus_n = (-d * (e - c) - sqrt(discriminator));
	float t_plus_n = (-d * (e - c) + sqrt(discriminator));
	float t_d = d * d;

	float t_minus = t_minus_n / t_d;
	float t_plus = t_plus_n / t_d;
	t = t_minus > 0 ? t_minus : t_plus;

	return t > 0;
}


Vector Sphere::getNormal( Vector point )
{
	Vector normal = point - center;
	return (normal.normalize());
}

AABB Sphere::GetBoundingBox() {
	Vector a_min = Vector(center.x-radius, center.y- radius, center.z - radius);
	Vector a_max = Vector(center.x + radius, center.y + radius, center.z + radius);
	return(AABB(a_min, a_max));
}

aaBox::aaBox(Vector& minPoint, Vector& maxPoint) //Axis aligned Box: another geometric object
{
	this->min = minPoint;
	this->max = maxPoint;
}

AABB aaBox::GetBoundingBox() {
	return(AABB(min, max));
}

bool aaBox::intercepts(Ray& ray, float& t)
{	
	double ox = ray.origin.x;
	double oy = ray.origin.y;
	double oz = ray.origin.z;
	double dx = ray.direction.x;
	double dy = ray.direction.y;
	double dz = ray.direction.z;
	double tx_min, ty_min, tz_min;
	double tx_max, ty_max, tz_max;

	double a = 1.0 / dx;
	if (a >= 0) {
		tx_min = (min.x - ox) * a;
		tx_max = (max.x - ox) * a;
	}
	else {
		tx_min = (max.x - ox) * a;
		tx_max = (min.x - ox) * a;
	}
	double b = 1.0 / dy;
	if (b >= 0) {
		ty_min = (min.y - oy) * b;
		ty_max = (max.y - oy) * b;
	}
	else {
		ty_min = (max.y - oy) * b;
		ty_max = (min.y - oy) * b;
	}
	double c = 1.0 / dz;
	if (c >= 0) {
		tz_min = (min.z - oz) * c;
		tz_max = (max.z - oz) * c;
	}
	else {
		tz_min = (max.z - oz) * c;
		tz_max = (min.z - oz) * c;
	}

	Vector face_in, face_out;
	float tE, tL; //entering and leaving t values Vector face_in, face_out; // normals
// find largest tE, entering t value
	if (tx_min > ty_min) {
		tE = tx_min;
		face_in = (a >= 0.0) ? Vector(-1, 0, 0) : Vector(1, 0, 0);
	}
	else {
		tE = ty_min;
		face_in = (b >= 0.0) ? Vector(0, -1, 0) : Vector(0, 1, 0);
	}
	if (tz_min > tE) { 
		tE = tz_min;
		face_in = (c >= 0.0) ? Vector(0, 0, -1) : Vector(0, 0, 1);
	}
	// find smallest tL, leaving t value
	if (tx_max < ty_max) {
		tL = tx_max;
		face_out = (a >= 0.0) ? Vector(1, 0, 0) : Vector(-1, 0, 0);
	}
	else {
		tL = ty_max;
		face_out = (b >= 0.0) ? Vector(0, 1, 0) : Vector(0, -1, 0);
	}
	if (tz_max < tL) { 
		tL = tz_max;
		face_out = (c >= 0.0) ? Vector(0, 0, 1) : Vector(0, 0, -1);
	}
	if (tE < tL && tL > 0) { // condition for a hit
		if (tE > 0) {
			t = tE; // ray hits outside surface
			Normal = face_in;
		}
		else {
			t = tL; // ray hits inside surface
			Normal = face_out;
		}
		return true;
	}
	return false;
}

Vector aaBox::getNormal(Vector point)
{
	return Normal;
}

Ray Physics::reflection(HitRecord hit, Vector& d) {
	Vector n = hit.frontFace ? hit.n : -hit.n;
	Vector r = d - n * (d * n) * 2;
	return Ray(hit.p + n * 0.0001, (r + sample_unit_sphere() * hit.object->GetMaterial()->GetRoughness()).normalize());
}

Ray Physics::refraction(HitRecord hit, Vector& d, float eta_i, float eta_t, float* reflection) {
	Vector n = hit.frontFace ? hit.n : -hit.n;
	Vector p_ = hit.frontFace ? hit.p  - hit.n * 0.0001 : hit.p + hit.n * 0.0001;
	float discriminator = 1 - ((1 - pow(d * n, 2)) * pow(eta_i, 2)) / pow(eta_t, 2);
	if (discriminator < 0) {
		*reflection = 1; // total reflection
		return Ray(p_, Vector(0,0,0));
	}
	else {
		Vector t_ = ((d - n * (d * n)) * eta_i) / eta_t - n * sqrt(discriminator);
		*reflection = Physics::reflectivity(hit, d, t_, eta_i, eta_t);
		return Ray(p_, t_);
	}
}

float Physics::reflectivity(HitRecord hit, Vector& d, Vector& t, float eta_i, float eta_t) {
	float r_0 = pow((eta_i - eta_t) / (eta_t + eta_i), 2);
	float cosTheta = hit.frontFace ? (-d * hit.n) : (-d * -hit.n);
	return r_0 + (1 - r_0) * pow(1 - cosTheta, 5);
}

Scene::Scene(Accelerator accelerator_, int decomposeLights_){
	accelerator = accelerator_;
	decomposeLights = decomposeLights_;
}

Scene::~Scene()
{
	/*for ( int i = 0; i < objects.size(); i++ )
	{
		delete objects[i];
	}
	objects.erase();
	*/
}

int Scene::getNumObjects()
{
	return objects.size();
}


void Scene::addObject(Object* o)
{
	objects.push_back(o);
}


Object* Scene::getObject(unsigned int index)
{
	if (index >= 0 && index < objects.size())
		return objects[index];
	return NULL;
}


int Scene::getNumLights()
{
	return lights.size();
}


void Scene::addLight(Light* l)
{
	lights.push_back(l);
}


Light* Scene::getLight(unsigned int index)
{
	if (index >= 0 && index < lights.size())
		return lights[index];
	return NULL;
}

Accelerator Scene::getAcceleration() {
	return accelerator;
}


void Scene::LoadSkybox(const char *sky_dir)
{
	char *filenames[6];
	char buffer[100];
	const char *maps[] = { "/right.jpg", "/left.jpg", "/top.jpg", "/bottom.jpg", "/front.jpg", "/back.jpg" };

	for (int i = 0; i < 6; i++) {
		strcpy_s(buffer, sizeof(buffer), sky_dir);
		strcat_s(buffer, sizeof(buffer), maps[i]);
		filenames[i] = (char *)malloc(sizeof(buffer));
		strcpy_s(filenames[i], sizeof(buffer), buffer);
	}
	
	ILuint ImageName;

	ilEnable(IL_ORIGIN_SET);
	ilOriginFunc(IL_ORIGIN_LOWER_LEFT);

	for (int i = 0; i < 6; i++) {
		ilGenImages(1, &ImageName);
		ilBindImage(ImageName);

		if (ilLoadImage(filenames[i]))  //Image loaded with lower left origin
			printf("Skybox face %d: Image sucessfully loaded.\n", i);
		else
			exit(0);

		ILint bpp = ilGetInteger(IL_IMAGE_BITS_PER_PIXEL);

		ILenum format = IL_RGB;
		printf("bpp=%d\n", bpp);
		if (bpp == 24)
			format = IL_RGB;
		else if (bpp == 32)
			format = IL_RGBA;

		ilConvertImage(format, IL_UNSIGNED_BYTE);

		int size = ilGetInteger(IL_IMAGE_SIZE_OF_DATA);
		skybox_img[i].img = (ILubyte *)malloc(size);
		ILubyte *bytes = ilGetData();
		memcpy(skybox_img[i].img, bytes, size);
		skybox_img[i].resX = ilGetInteger(IL_IMAGE_WIDTH);
		skybox_img[i].resY = ilGetInteger(IL_IMAGE_HEIGHT);
		format == IL_RGB ? skybox_img[i].BPP = 3 : skybox_img[i].BPP = 4;
		ilDeleteImages(1, &ImageName);
	}
	ilDisable(IL_ORIGIN_SET);
}

Color Scene::GetSkyboxColor(Ray& r) {
	float t_intersec;
	Vector cubemap_coords; //To index the skybox

	float ma;
	CubeMap img_side;
	float sc, tc, s, t;
	unsigned int xp, yp, width, height, bytesperpixel;

	//skybox indexed by the ray direction
	cubemap_coords = r.direction;


	if (fabs(cubemap_coords.x) > fabs(cubemap_coords.y)) {
		ma = fabs(cubemap_coords.x);
		cubemap_coords.x >= 0 ? img_side = LEFT : img_side = RIGHT;    //left cubemap at X = +1 and right at X = -1
	}
	else {
		ma = fabs(cubemap_coords.y);
		cubemap_coords.y >= 0 ? img_side = TOP : img_side = BOTTOM; //top cubemap at Y = +1 and bottom at Y = -1
	}

	if (fabs(cubemap_coords.z) > ma) {
		ma = fabs(cubemap_coords.z);
		cubemap_coords.z >= 0 ? img_side = FRONT : img_side = BACK;   //front cubemap at Z = +1 and back at Z = -1
	}

	switch (img_side) {

	case 0:  //right
		sc = -cubemap_coords.z;
		tc = cubemap_coords.y;
		break;

	case 1:  //left
		sc = cubemap_coords.z;
		tc = cubemap_coords.y;
		break;

	case 2:  //top
		sc = -cubemap_coords.x;
		tc = -cubemap_coords.z;
		break;

	case 3: //bottom
		sc = -cubemap_coords.x;
		tc = cubemap_coords.z;
		break;

	case 4:  //front
		sc = -cubemap_coords.x;
		tc = cubemap_coords.y;
		break;

	case 5: //back
		sc = cubemap_coords.x;
		tc = cubemap_coords.y;
		break;
	}

	double invMa = 1 / ma;
	s = (sc * invMa + 1) / 2;
	t = (tc * invMa + 1) / 2;

	width = skybox_img[img_side].resX;
	height = skybox_img[img_side].resY;
	bytesperpixel = skybox_img[img_side].BPP;

	xp = int((width - 1) * s);
	xp < 0 ? 0 : (xp > (width - 1) ? width - 1 : xp);
	yp = int((height - 1) * t);
	yp < 0 ? 0 : (yp > (height - 1) ? height - 1 : yp);

	float red = u8tofloat(skybox_img[img_side].img[(yp*width + xp) * bytesperpixel]);
	float green = u8tofloat(skybox_img[img_side].img[(yp*width + xp) * bytesperpixel + 1]);
	float blue = u8tofloat(skybox_img[img_side].img[(yp*width + xp) * bytesperpixel + 2]);

	return(Color(red, green, blue));
}




////////////////////////////////////////////////////////////////////////////////
// P3F file parsing methods.
//
void next_token(ifstream& file, char *token, const char *name)
{
  file >> token;
  if (strcmp(token, name))
    cerr << "'" << name << "' expected.\n";
}

bool Scene::load_p3f(const char *name)
{
  const	int	lineSize = 1024;
  string	cmd;
  char		token	[256];
  ifstream	file(name, ios::in);
  Material *	material;

  material = NULL;

  if (file >> cmd)
  {
    while (true)
    {
      
	  if (cmd == "f")   //Material
      {
	    double Kd, Ks, Shine, T, ior, roughness = 0;
	    Color cd, cs;
		string s;
		getline(file, s);
		istringstream iss(s);

		iss >> cd >> Kd >> cs >> Ks >> Shine >> T >> ior >> roughness;

		material = new Material(cd, Kd, cs, Ks, Shine, T, ior, roughness);
      }

      else if (cmd == "s")    //Sphere
      {
	     Vector center;
    	 float radius;
         Sphere* sphere;

	    file >> center >> radius;
        sphere = new Sphere(center,radius);
	    if (material) sphere->SetMaterial(material);
        this->addObject( (Object*) sphere);
      }

	  else if (cmd == "box")    //axis aligned box
	  {
		  Vector minpoint, maxpoint;
		  aaBox	*box;

		  file >> minpoint >> maxpoint;
		  box = new aaBox(minpoint, maxpoint);
		  if (material) box->SetMaterial(material);
		  this->addObject((Object*)box);
	  }
	  else if (cmd == "p")  // Polygon: just accepts triangles for now
      {
		  Vector P0, P1, P2;
		  Triangle* triangle;
		  unsigned total_vertices;
		  
		  file >> total_vertices;
		  if (total_vertices == 3)
		  {
			  file >> P0 >> P1 >> P2;
			  triangle = new Triangle(P0, P1, P2);
			  if (material) triangle->SetMaterial(material);
			  this->addObject( (Object*) triangle);
		  }
		  else
		  {
			  cerr << "Unsupported number of vertices.\n";
			  break;
		  }
      }
      
	  else if (cmd == "mesh") {
		  unsigned total_vertices, total_faces;
		  unsigned P0, P1, P2;
		  Triangle* triangle;
		  Vector* verticesArray, vertex;

		  file >> total_vertices >> total_faces;
		  verticesArray = (Vector*)malloc(total_vertices * sizeof(Vector));
		  for (int i = 0; i < total_vertices; i++) {
			  file >> vertex;
			  verticesArray[i] = vertex;
		  }
		  for (int i = 0; i < total_faces; i++) {
			  file >> P0 >> P1 >> P2;
			  triangle = new Triangle(verticesArray[P0 - 1], verticesArray[P1 - 1], verticesArray[P2 - 1]); //vertex index start at 1
			  if (material) triangle->SetMaterial(material);
			  this->addObject((Object*)triangle);
		  }

	  }

	  else if (cmd == "pl")  // General Plane
	  {
          Vector P0, P1, P2;
		  Plane* plane;

          file >> P0 >> P1 >> P2;
          plane = new Plane(P0, P1, P2);
	      if (material) plane->SetMaterial(material);
          this->addObject( (Object*) plane);
	  }

      else if (cmd == "l")  // Need to check light color since by default is white
      {
	    Vector pos;
        Color color;

	    file >> pos >> color;
	    
	      this->addLight(new Light(pos, color));
	    
      }
	  else if (cmd == "al")
	  {
		  Vector pos;
		  Vector at;
		  float width, height;
		  Color color;

		  file >> pos >> at >> width >> height >> color;

		  AreaLight* light = new AreaLight(pos, color, at, width, height);

		  if (decomposeLights > 1) {
			  vector<Light*> lights = light->decompose(decomposeLights);
			  for (Light* pointLight : lights) {
				  this->addLight(pointLight);
			  }
		  }
		  else {
			  this->addLight(light);
		  }
	  }
      else if (cmd == "v")
      {
	    Vector up, from, at;
	    float fov, hither;
	    int xres, yres;
        Camera* camera;
		float focal_ratio; //ratio beteween the focal distance and the viewplane distance
		float aperture_ratio; // number of times to be multiplied by the size of a pixel

	    next_token (file, token, "from");
	    file >> from;

	    next_token (file, token, "at");
	    file >> at;

	    next_token (file, token, "up");
	    file >> up;

	    next_token (file, token, "angle");
	    file >> fov;

	    next_token (file, token, "hither");
	    file >> hither;

	    next_token (file, token, "resolution");
	    file >> xres >> yres;

		next_token(file, token, "aperture");
		file >> aperture_ratio;

		next_token(file, token, "focal");
		file >> focal_ratio;
	    // Create Camera
		camera = new Camera( from, at, up, fov, hither, 100.0*hither, xres, yres, aperture_ratio, focal_ratio);
        this->SetCamera(camera);
      }

      else if (cmd == "bclr")   //Background color
      {
		Color bgcolor;
		file >> bgcolor;
		this->SetBackgroundColor(bgcolor);
	  }
	
	  else if (cmd == "env")
	  {
		  file >> token;
		  
		  this->LoadSkybox(token);
		  this->SetSkyBoxFlg(true);
	  }
      else if (cmd[0] == '#')
      {
	    file.ignore (lineSize, '\n');
      }
      else
      {
	    cerr << "unknown command '" << cmd << "'.\n";
	    break;
      }
      if (!(file >> cmd))
        break;
    }
  }

  file.close();
  return true;
};

void Scene::create_random_scene() {
	Camera* camera;
	Material* material;
	Sphere* sphere;

	set_rand_seed(time(NULL) * time(NULL) * time(NULL));
	material = NULL;
	this->SetSkyBoxFlg(false);  //init with no skybox

	this->SetBackgroundColor(Color(0.5, 0.7, 1.0));
	//this->LoadSkybox("skybox");
	//this->SetSkyBoxFlg(true);
	
	camera = new Camera(Vector(13.0, 2.0, 3.0), Vector(0.0, 0.0, 0), Vector(0.0, 1.0, 0.0), 45.0, 0.01, 10000.0, 800, 600, 0, 1.5f);
	this->SetCamera(camera);

	this->addLight(new Light(Vector(7, 10, -5), Color(1.0, 1.0, 1.0)));
	this->addLight(new Light(Vector(-7, 10, -5), Color(1.0, 1.0, 1.0)));
	this->addLight(new Light(Vector(0, 10, 7), Color(1.0, 1.0, 1.0)));

	material = new Material(Color(0.5, 0.5, 0.5), 1.0, Color(0.0, 0.0, 0.0), 0.0, 10, 0, 1);


	sphere = new Sphere(Vector(0.0, -1000, 0.0), 1000.0);
	if (material) sphere->SetMaterial(material);
	this->addObject((Object*)sphere);

	for (int a = -5; a < 5; a++)
		for (int b = -5; b < 5; b++) {

			double choose_mat = rand_double();

			Vector center = Vector(a + 0.9 * rand_double(), 0.2, b + 0.9 * rand_double());

			if ((center - Vector(4.0, 0.2, 0.0)).length() > 0.9) {
				if (choose_mat < 0.4) {  //diffuse
					material = new Material(Color(rand_double(), rand_double(), rand_double()), 1.0, Color(0.0, 0.0, 0.0), 0.0, 10, 0, 1);
					sphere = new Sphere(center, 0.2);
					if (material) sphere->SetMaterial(material);
					this->addObject((Object*)sphere);
				}
				else if (choose_mat < 0.9) {   //metal
					material = new Material(Color(0.0, 0.0, 0.0), 0.0, Color(rand_double(0.5, 1), rand_double(0.5, 1), rand_double(0.5, 1)), 1.0, 220, 0, 1);
					sphere = new Sphere(center, 0.2);
					if (material) sphere->SetMaterial(material);
					this->addObject((Object*)sphere);
				}
				else {   //glass 
					material = new Material(Color(0.0, 0.0, 0.0), 0.0, Color(1.0, 1.0, 1.0), 0.7, 20, 1, 1.5);
					sphere = new Sphere(center, 0.2);
					if (material) sphere->SetMaterial(material);
					this->addObject((Object*)sphere);
				}

			}

		}

	material = new Material(Color(0.0, 0.0, 0.0), 0.0, Color(1.0, 1.0, 1.0), 0.7, 20, 1, 1.5);
	sphere = new Sphere(Vector(0.0, 1.0, 0.0), 1.0);
	if (material) sphere->SetMaterial(material);
	this->addObject((Object*)sphere);

	material = new Material(Color(0.4, 0.2, 0.1), 0.9, Color(1.0, 1.0, 1.0), 0.1, 10, 0, 1.0);
	sphere = new Sphere(Vector(-4.0, 1.0, 0.0), 1.0);
	if (material) sphere->SetMaterial(material);
	this->addObject((Object*)sphere);

	material = new Material(Color(0.4, 0.2, 0.1), 0.0, Color(0.7, 0.6, 0.5), 1.0, 220, 0, 1.0);
	sphere = new Sphere(Vector(4.0, 1.0, 0.0), 1.0);
	if (material) sphere->SetMaterial(material);
	this->addObject((Object*)sphere);
}

void Scene::build() {
	vector<Object*> objs;
	int num_objects = getNumObjects();
	for (int o = 0; o < num_objects; o++) {
		objs.push_back(getObject(o));
	}

	if (accelerator == GRID_ACC) {
		grid = new Grid();
		grid->Build(objs);
	}
	else if (accelerator == BVH_ACC) {
		bvh = new BVH();
		bvh->Build(objs);
	}
	printf("Scene built.\n\n");
}

bool Scene::traverseScene(Ray& ray, Object** object, Vector& hitpoint) {
	if (accelerator == GRID_ACC) {
		return grid->Traverse(ray, object, hitpoint);
	}
	else if (accelerator == BVH_ACC) {
		return bvh->Traverse(ray, object, hitpoint);
	}
	return false;
}

bool Scene::traverseSceneShadow(Ray& ray) {
	if (accelerator == GRID_ACC) {
		return grid->Traverse(ray);
	}
	else if (accelerator == BVH_ACC) {
		return bvh->Traverse(ray);
	}
	return false;
}

bool Scene::traverseGrid(Ray& ray, Object** object, Vector& hitpoint) {
	return grid->Traverse(ray, object, hitpoint);
}

bool Scene::traverseShadowGrid(Ray& ray) {
	return grid->Traverse(ray);
}