inline float poly_6(float r, float h) {
	return 
		(1.f - clamp(floor(r / h), 0.f, 1.f)) *
		(315.f / (64.f * M_PI * pown(h, 9))) * pown((pown(h, 2) - pown(r, 2)), 3);
}

inline float3 poly_6_gradient(float3 r, float h) {
	return
		(1.f - clamp(floor(convert_float(length(r)) / h), 0.f, 1.f)) *
		convert_float(-945.f / (32.f * M_PI * pown(h, 9))) * 
		r * 
		pown((pown(h, 2) - pown(length(r), 2)), 2);
}

inline float poly_6_laplacian(float r, float h) {
	return
		(1.f - clamp(floor(convert_float(length(r)) / h), 0.f, 1.f)) *
		-945.f / (32.f * M_PI * pown(h, 9)) * 
		(pown(h, 2) - pown(r, 2)) * 
		(3.f * pown(h, 2) - 7.f * pown(r, 2));
}

#define EPSILON 0.0000001f

inline float3 spiky_gradient(float3 r, float h) {
	if(length(r) - EPSILON < 0.f && 0.f < length(r) + EPSILON) {
		return (-45.f / convert_float((M_PI * pown(h, 6))));
	}
	return
	 	(1.f - clamp(floor(convert_float(length(r)) / h), 0.f, 1.f)) *
		(-45.f / convert_float((M_PI * pown(h, 6)))) * 
		(r / convert_float(length(r))) * 
		pown(h - convert_float(length(r)), 2);
}

inline float viscosity_laplacian(float r, float h) {
	return 
		(1.f - clamp(floor(r / h), 0.f, 1.f)) *
		(45.f / (M_PI * pown(h, 6))) * (h - r);
}