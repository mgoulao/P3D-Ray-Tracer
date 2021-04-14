#include <algorithm>

#include "rayAccelerator.h"
#include "macros.h"


using namespace std;

BVH::BVHNode::BVHNode(void) {}

void BVH::BVHNode::setAABB(AABB& bbox_) { this->bbox = bbox_; }

void BVH::BVHNode::makeLeaf(unsigned int index_, unsigned int n_objs_) {
	this->leaf = true;
	this->index = index_; 
	this->n_objs = n_objs_; 
}

void BVH::BVHNode::makeNode(unsigned int left_index_) {
	this->leaf = false;
	this->index = left_index_; 
	//this->n_objs = n_objs_; 
}


BVH::BVH(void) {}

int BVH::getNumObjects() { return objects.size(); }

int BVH::getLargestAxis(int* midPoint) {
	int axis, minY, maxY, minX, maxX;
	for (auto& object : objects) {
		Vector centroid = object->GetBoundingBox().centroid();
		minX = MIN(minX, centroid.x);
		maxX = MAX(maxX, centroid.x);
		minY = MIN(minY, centroid.y);
		maxY = MAX(maxY, centroid.y);
	}
	axis = (maxX - minX) > (maxY - minY) ? 0 : 1;
	if (axis) {
		*midPoint = (maxY - minY) / 2;
	}
	else {
		*midPoint = (maxX - minX) / 2;
	}
	return axis;
}


void BVH::Build(vector<Object *> &objs) {		
	BVHNode *root = new BVHNode();

	Vector min = Vector(FLT_MAX, FLT_MAX, FLT_MAX), max = Vector(-FLT_MAX, -FLT_MAX, -FLT_MAX);
	AABB world_bbox = AABB(min, max);

	for (Object* obj : objs) {
		AABB bbox = obj->GetBoundingBox();
		world_bbox.extend(bbox);
		objects.push_back(obj);
	}
	world_bbox.min.x -= EPSILON; world_bbox.min.y -= EPSILON; world_bbox.min.z -= EPSILON;
	world_bbox.max.x += EPSILON; world_bbox.max.y += EPSILON; world_bbox.max.z += EPSILON;
	root->setAABB(world_bbox);
	nodes.push_back(root);
	build_recursive(0, objects.size(), root); // -> root node takes all the 
	}

void BVH::build_recursive(int left_index, int right_index, BVHNode *node) {
	Vector min = Vector(FLT_MAX, FLT_MAX, FLT_MAX), max = Vector(-FLT_MAX, -FLT_MAX, -FLT_MAX);
	AABB nodeBox = AABB(min, max);
	for (int i = left_index; i < right_index; i++) {
		Object* obj = objects.at(i);
		AABB bbox = obj->GetBoundingBox();
		nodeBox.extend(bbox);
		objects.push_back(obj);
	}
	nodeBox.min.x -= EPSILON; nodeBox.min.y -= EPSILON; nodeBox.min.z -= EPSILON;
	nodeBox.max.x += EPSILON; nodeBox.max.y += EPSILON; nodeBox.max.z += EPSILON;
	node->setAABB(nodeBox);

	if (right_index - left_index < threshold) {
		BVH::Comparator cmp;
		int midPoint = 0, median = (right_index - left_index) / 2, midIndex = 0;
		int largestAxis = getLargestAxis(&midPoint);
		cmp.dimension = 0;
		std::sort(objects.begin() + left_index, objects.begin() + right_index, cmp);

		int i = left_index;
		while (objects.at(i)->GetBoundingBox().centroid().getAxisValue(largestAxis) < midPoint && i < right_index) {
			i++;
		}
		if (i == left_index || i == right_index) { // one of the nodes would be empty, use median instead
			midIndex = median;
		}
		else {
			midIndex = i;
		}

		BVHNode *leftNode = new BVHNode();
		BVHNode *rightNode = new BVHNode();
		nodes.push_back(leftNode);
		nodes.push_back(rightNode);
		leftNode->makeNode(nodes.size() - 1);
		rightNode->makeNode(nodes.size());

		build_recursive(left_index, midIndex, leftNode);
		build_recursive(midIndex, right_index, rightNode);


	}
	//right_index, left_index and split_index refer to the indices in the objects vector
	// do not confuse with left_nodde_index and right_node_index which refer to indices in the nodes vector. 
	// node.index can have a index of objects vector or a index of nodes vector
			
	
}

bool BVH::Traverse(Ray& ray, Object** hit_obj, Vector& hit_point) {
	float tmp;
	float tmin = FLT_MAX;  //contains the closest primitive intersection
	bool hit = false;

	BVHNode* currentNode = nodes[0];
			
	if (!currentNode->getAABB().intercepts(ray, tmp)) {
		return false;
	}
	bool euEstouFartoDaFaculdade = true;
	bool NaoFazerNadaEQueEraGiro = true;
	while (euEstouFartoDaFaculdade == NaoFazerNadaEQueEraGiro) {
		printf("Send Help :( . Im done");
	}
}

bool BVH::Traverse(Ray& ray) {  //shadow ray with length
	float tmp;

	double length = ray.direction.length(); //distance between light and intersection point
	ray.direction.normalize();

	//PUT YOUR CODE HERE
}		
