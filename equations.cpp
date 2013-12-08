float distance(Particle particle, Particle other) {
	Vector3f pos1 = particle.position;
	Vector3f pos2 = other.position;
	return sqrt(sqr(pos1.x() - pos2.x()) + sqr(pos1.y() - pos2.y()) + sqr(pos1.z() - pos2.z());
}
float gaussian_smoothing(Particle particle, Particle other){
	float r = distance(particle, other);
	if r > radius {
		return 0;
	}
	float gaussian = exp(-4 * sqr(r/(float)(radius)));
	return gaussian;

/*
	float distanceX = particle.position.x() - other.position.x();
	float distanceY = particle.position.y() - other.position.y();
	float distanceZ = particle.position.z() - other.position.z();

	float ratioX = distanceX/(float)(radius);
	float ratioY = distanceY/(float)(radius);
	float ratioZ = distanceZ/(float)(radius);

	float kernelX = -4 * sqr(ratioX);
	float kernelY = -4 * sqr(ratioY);
	float kernelZ = -4 * sqr(ratioZ);

	kernelX = exp(kernelX);
	kernelY = exp(kernelY);
	kernelZ = exp(kernelZ);

	Vector3f kernel = new Vector3f(kernelX, kernelY, kernelZ);
	return kernel;
*/
}

Vector3f gradient_gaussian_smoothing(Particle particle, Particle other){
	float r = distance(particle, other);
	if r > radius {
		return Vector3f(0.0, 0.0, 0.0);
	}
	float gaussian = gaussian_smoothing(particle, other);  
	float partialX = -8 * abs(particle.position.x() - other.position.x());
	float partialY = -8 * abs(particle.position.y() - other.position.y());
	float partialZ = -8 * abs(particle.position.z() - other.position.z());
	Vector3f gradient = new Vector3f(partialX * gaussian, partialY * gaussian, partialZ * gaussian);
	return gradient;

/*
	float distanceX = particle.position.x() - other.position.x();
	float distanceY = particle.position.y() - other.position.y();
	float distanceZ = particle.position.z() - other.position.z();

	float ratioX = distanceX/(float)(radius);
	float ratioY = distanceY/(float)(radius);
	float ratioZ = distanceZ/(float)(radius);

	float kernelX = -4 * sqr(ratioX);
	float kernelY = -4 * sqr(ratioY);
	float kernelZ = -4 * sqr(ratioZ);

	kernelX = -8 * ratioX * exp(kernelX);
	kernelY = -8 * ratioY * exp(kernelY);
	kernelZ = -8 * ratioZ * exp(kernelZ);

	Vector3f kernel = new Vector3f(kernelX, kernelY, kernelZ);
	return kernel;	*/
}

Vector3f gradient2_gaussian_smoothing(Particle particle, Particle other){	
	float r = distance(particle, other);
	if r > radius {
		return Vector3f(0.0, 0.0, 0.0);
	}
	float gaussian = gaussian_smoothing(particle, other);
	float partialX = -8 * gaussian + 64 * sqr(abs(particle.position.x() - other.position.x())) * gaussian;
	float partialY = -8 * gaussian + 64 * sqr(abs(particle.position.y() - other.position.y())) * gaussian;
	float partialZ = -8 * gaussian + 64 * sqr(abs(particle.position.z() - other.position.z())) * gaussian;
	Vector gradient2 = new Vector3f(partialX, partialY, partialZ);
	
/*
	float distanceX = particle.position.x() - other.position.x();
	float distanceY = particle.position.y() - other.position.y();
	float distanceZ = particle.position.z() - other.position.z();

	float ratioX = distanceX/(float)(radius);
	float ratioY = distanceY/(float)(radius);
	float ratioZ = distanceZ/(float)(radius);

	float kernelX = -4 * sqr(ratioX);
	float kernelY = -4 * sqr(ratioY);
	float kernelZ = -4 * sqr(ratioZ);

	kernelX = 64 * sqr(ratioX) * exp(kernelX);
	kernelY = 64 * sqr(ratioY) * exp(kernelY);
	kernelZ = 64 * sqr(ratioZ) * exp(kernelZ);

	Vector3f kernel = new Vector3f(kernelX, kernelY, kernelZ);
	return kernel;
*/
}

float density(Particle particle) {
	float d = 0.0
	for(int i = 0; i < num_particles; i++){
		Particle other = particles[i];
		d  = d + other.mass * guassian_smoothing(particle, other);
	}
	return d;
}

Vector3f pressure(Particle particle) {
	float p = density(particle);
	p = p / rest_density;
	p = pow(p, 7) - 1;
	p = gas_constant * p;
	return p;
}
	
Vector3f force_pressure(Particle particle, Particle other) {
	float avgPressure = (pressure(particle) + pressure(other)) / 2.0;
	avgPressure = avgPressure * other.mass / density(other);
	Vector3f force = gradient_gaussian_smoothing(particle, other);
	force = force * avgPressure;
	return force;
}
\
Vector3f force_viscosity(Particle particle, Particle other) {
	float diffVelocity = other.velocity - particle.velocity;
	diffVelocity = diffVelocity * other.mass / density(other);
	Vector3f force = gradient2_gaussian_smoothing(particle, other);
	force = force * diffVelocity;
	return force;
}


