#include <algorithm>

#include "rayAccelerator.h"
#include "macros.h"


using namespace std;

bool BVH::Comparator::operator() (Object* a, Object* b) {
	float ca = a->GetBoundingBox().centroid().getAxisValue(dimension);
	float cb = b->GetBoundingBox().centroid().getAxisValue(dimension);
	return ca < cb;
}

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
}


BVH::BVH(void) {}

int BVH::getLargestAxis(int* midPoint) {
	int axis, minY = INT_MAX, maxY = -INT_MAX, minX = INT_MAX, maxX = -INT_MAX, 
		minZ = INT_MAX, maxZ = -INT_MAX;
	for (auto& object : objects) {
		Vector centroid = object->GetBoundingBox().centroid();
		minX = MIN(minX, centroid.x);
		maxX = MAX(maxX, centroid.x);
		minY = MIN(minY, centroid.y);
		maxY = MAX(maxY, centroid.y);
		minZ = MIN(minZ, centroid.z);
		maxZ = MAX(maxZ, centroid.z);
	}
	int dX = maxX - minX, dY = maxY - minY, dZ = maxZ - minZ;
	axis = (dX > dY ? (dX > dZ ? 0 : 2) : (dY > dZ ? 1 : 2));
	if (axis == 0) {
		*midPoint = minX + (dX / 2);
	}
	else if (axis == 1) {
		*midPoint = minY + (dY / 2);
	}
	else {
		*midPoint = minZ + (dZ / 2);
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
		addObject(obj);
	}

	world_bbox.min.x -= EPSILON; world_bbox.min.y -= EPSILON; world_bbox.min.z -= EPSILON;
	world_bbox.max.x += EPSILON; world_bbox.max.y += EPSILON; world_bbox.max.z += EPSILON;
	root->setAABB(world_bbox);
	nodes.push_back(root);
	root->makeNode(nodes.size());
	build_recursive(0, objects.size(), root); // -> root node takes all the 
}

void BVH::setNodeAABB(BVHNode* node, int left_index, int right_index) {
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
}

void BVH::build_recursive(int left_index, int right_index, BVHNode *node) {

	if (right_index - left_index > threshold) {
		BVH::Comparator cmp;
		int midPoint = 0, median = left_index + (right_index - left_index) / 2, midIndex = 0;
		int largestAxis = getLargestAxis(&midPoint);
		int i = left_index;
		
		cmp.dimension = largestAxis;
		std::sort(objects.begin() + left_index, objects.begin() + right_index, cmp);
		node->makeNode(nodes.size());

		while (objects.at(i)->GetBoundingBox().centroid().getAxisValue(largestAxis) < midPoint && i < right_index-1) {
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

		setNodeAABB(leftNode, left_index, midIndex);
		setNodeAABB(rightNode, midIndex, right_index);

		build_recursive(left_index, midIndex, leftNode);
		build_recursive(midIndex, right_index, rightNode);
	}
	else { // leaf
		node->makeLeaf(left_index, right_index - left_index);
		setNodeAABB(node, left_index, right_index);
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
		BVHNode* closestNode = NULL;
		if (!currentNode->isLeaf()) {
			BVHNode* leftNode = nodes.at(currentNode->getLeftNodeIndex());
			BVHNode* rightNode = nodes.at(currentNode->getRightNodeIndex());
			float leftT = 0, rightT = 0;
			bool leftHit = leftNode->getAABB().intercepts(ray, leftT);
			bool rightHit = rightNode->getAABB().intercepts(ray, rightT);

			if (leftHit && rightHit) {
				closestNode = leftT <= rightT ? leftNode : rightNode;
				BVHNode* farthestNode = leftT <= rightT ? rightNode : leftNode;
				float farthestT = MAX(leftT, rightT);
				hit_stack.push(StackItem(farthestNode, farthestT));

				currentNode = closestNode;
			}
			else if (leftHit || rightHit) {
				closestNode = leftHit ? leftNode : rightNode;
				currentNode = closestNode;
			}
		}
		else {
			int firstObjectIndex = currentNode->getIndex();
			for (int i = 0; i < currentNode->getNObjs(); i++) {
				float t = 0;
				Object* obj = objects.at(firstObjectIndex + i);
				bool objHit = obj->intercepts(ray, t);
				hit = hit || objHit;
				if (objHit && t < tmin) {
					*hit_obj = obj;
					hit_point = ray.origin + ray.direction * t;
					tmin = t;
				}
			}
		}

		if (closestNode == NULL) {
			if (hit_stack.size()) {
				bool empty = false;
				StackItem stackItem = hit_stack.top();
				while (stackItem.t > tmin) {
					hit_stack.pop();
					if (!hit_stack.size()) {
						empty = true;
						break;
					}
					stackItem = hit_stack.top();
				}
				if (empty) {
					break;
				}
				currentNode = stackItem.ptr;
				hit_stack.pop();
			} else {
				break;
			}
		}
	}
	return hit;
}


bool BVH::Traverse(Ray& ray) {  //shadow ray with length
	float tmp;
	float length = ray.direction.length(); //distance between light and intersection point
	ray.direction.normalize();

	BVHNode* currentNode = nodes[0];

	if (!currentNode->getAABB().intercepts(ray, tmp)) {
		return false;
	}

	bool euEstouFartoDaFaculdade = true;
	bool NaoFazerNadaEQueEraGiro = true;
	while (euEstouFartoDaFaculdade == NaoFazerNadaEQueEraGiro) {
		if (!currentNode->isLeaf()) {
			BVHNode* leftNode = nodes.at(currentNode->getLeftNodeIndex());
			BVHNode* rightNode = nodes.at(currentNode->getRightNodeIndex());
			bool leftHit = leftNode->getAABB().intercepts(ray, tmp);
			bool rightHit = rightNode->getAABB().intercepts(ray, tmp);

			if (leftHit && rightHit) { // hits both
				hit_stack.push(StackItem(rightNode, 0));
				currentNode = leftNode;
				continue;
			}
			else if (leftHit || rightHit) {
				currentNode = leftHit ? leftNode : rightNode;
				continue;
			}
		}
		else {
			int firstObjectIndex = currentNode->getIndex();
			for (int i = 0; i < currentNode->getNObjs(); i++) {
				float t = 0;
				Object* obj = objects.at(firstObjectIndex + i);
				bool objHit = obj->intercepts(ray, t);
				if (objHit && t < length) {
					return true;
				}
			}
		}

		if (hit_stack.size()) {
			StackItem stackItem = hit_stack.top();
			currentNode = stackItem.ptr;
			hit_stack.pop();
		}
		else {
			break;
		}
	}
	return false;
}		
