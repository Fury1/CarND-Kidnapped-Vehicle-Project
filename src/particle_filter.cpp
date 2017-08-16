/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <algorithm>
#include <numeric>
#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	// Total particle count to represent vehicle location and uncertainties space
	num_particles = 200;

	// Set an initial location distribution based upon GPS coordinates and measurement noise
	normal_distribution<double> normal_x(x, std[0]);
	normal_distribution<double> normal_y(y, std[1]);
	normal_distribution<double> normal_theta(theta, std[2]); // aka yaw

	// Create the particles from the distribution randomly and add them to the particles vector
	// Add the corresponding initialized weight to the weight vector
	for (int i = 0; i < num_particles; i++) {
		Particle particle;
		particle.id = i;
		particle.x = normal_x(gen);
		particle.y = normal_y(gen);
		particle.theta = normal_theta(gen);
		particle.weight = 1.0; // initialize all weights to 1

		// Add to respective vector
		particles.push_back(particle);
		weights.push_back(1.0); // initialize all weights to 1
	}

	// Set initialization to true now that we are initialized to rough GPS coordinates
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	// Move every particle randomly based on a bicycle motion model to precurse the next weight update
	for (int i = 0; i < num_particles; i++) {
		double new_x;
		double new_y;
		double new_theta;

		/*
		* Bicycle motion model random movement of particles
		* Check for divide by 0 aka a straight path
		*/
		if (yaw_rate == 0) {
			new_x = particles[i].x + velocity * delta_t * cos(particles[i].theta);
			new_y = particles[i].y + velocity * delta_t * sin(particles[i].theta);
			new_theta = particles[i].theta;
		}
		else {
			new_x = particles[i].x + (velocity / yaw_rate) * (sin(particles[i].theta + (yaw_rate * delta_t)) - sin(particles[i].theta));
			new_y = particles[i].y + (velocity / yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + (yaw_rate * delta_t)));
			new_theta = particles[i].theta + (yaw_rate * delta_t);
		}

		// Create a gaussian distribution from the new particle position
		normal_distribution<double> normal_x(new_x, std_pos[0]);
		normal_distribution<double> normal_y(new_y, std_pos[1]);
		normal_distribution<double> normal_theta(new_theta, std_pos[2]); // aka yaw

		// Update the particles based on the motion model values
		// Randomly select these values from the gaussian distribution of new_x and new_y
		particles[i].x = normal_x(gen);
		particles[i].y = normal_y(gen);
		particles[i].theta = normal_theta(gen); // aka yaw
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	/*
	* Transform all coordinates into the same coordinate system
	* Observation coordinates are transformed into the global coordinate space/map space
	*/

	// Calculate normalization term for weight update
	const double gauss_norm = 1 / (2 * M_PI * std_landmark[0] * std_landmark[1]);

	// Transform each particle and update its weights respectively
	for (int i = 0; i < num_particles; i++) {
		// Create a vector for observation weights used to create the final particle weight
		vector<double> observation_weights;
		// Verified transformed observation landmark associations
		vector<int> associations; 	// landmark id
		vector<double> sense_x; 	// landmark x
		vector<double> sense_y; 	// landmark y

		// Map the vehicle observations to the map coordinate space with respect to the particle
		for (int j = 0; j < observations.size(); j++) {

			double cos_p = cos(particles[i].theta); // cosine of theta
			double sin_p = sin(particles[i].theta); // sine of theta

			LandmarkObs transformed_observation = {
				// Sync the observation id with the transformed version being made to stay organized (just in case its needed)
				observations[i].id,
				// Tranform the x, y value for the current particle to map (global) space
				// Rotation and shift of coordinates
				particles[i].x + (observations[j].x * cos_p - observations[j].y * sin_p),
				particles[i].y + (observations[j].x * sin_p + observations[j].y * cos_p)
			};

			// Find the closest landmark to the current particle observation (nearest neighboor), calculate the weight
			double nearest_landmark_range = numeric_limits<double>::max();	 // range to closest landmark found, start with max value
			LandmarkObs nearest_landmark; 									// holds the closest landmark observation information

			// Check the transformed observation against the landmarks and update the nearest landmark and range (see above)
			for (int k = 0; k < map_landmarks.landmark_list.size(); k++) {

				// Abbreviate the landmark for convenience
				Map::single_landmark_s landmark = map_landmarks.landmark_list[k];

				// Find the distance of the observation to the current landmark
				double landmark_range = dist(transformed_observation.x,
											 transformed_observation.y,
											 landmark.x_f,
											 landmark.y_f);
				/*
				* If the range of the observation to a landmark is lower then the sensor range AND less then the current
				* nearest landmark, update nearest landmark and its range
				*/
				if (landmark_range < nearest_landmark_range && landmark_range < sensor_range) {
					nearest_landmark_range = landmark_range;
					nearest_landmark.id = landmark.id_i;
					nearest_landmark.x = landmark.x_f;
					nearest_landmark.y = landmark.y_f;
				}
				else {
					// If distance to landmark is over sensor range or greater then previous, move on to the next landmark
					continue;
				}
			}

			// Save the best observation landmark associations to update the particle later
			associations.push_back(nearest_landmark.id);
			sense_x.push_back(nearest_landmark.x);
			sense_y.push_back(nearest_landmark.y);

			// With the nearest landmark to the observation, calculate the weight and add it to the particle weight vector
			// Calculate exponent
			double exponent = (pow(transformed_observation.x - nearest_landmark.x, 2)) /
							  (2 * pow(std_landmark[0], 2)) +
							  (pow(transformed_observation.y - nearest_landmark.y, 2)) /
							  (2 * pow(2 * std_landmark[1], 2));

			// Calculate the weight for this observation
			double weight = gauss_norm * exp(-exponent);

			// Add the weight to weights vector for final particle weight calculation later
			observation_weights.push_back(weight);
		}

	// Set all of the landmark associations to the particle, this is for visualization purposes in the simulator only
	particles[i] = SetAssociations(particles[i], associations, sense_x, sense_y);

	// Calculate the particles final weight
	for (int l = 0; l < observation_weights.size(); l++) {
		particles[i].weight *= observation_weights[l];
	}

	// Update the particle filters weight vector
	weights[i] = particles[i].weight;
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	// Create a normalized weighted index based on the particle weights, higher weights have more index occurences
    discrete_distribution<int> distribution(weights.begin(), weights.end());
	vector<Particle> resampled_particles; // new vector of particles

	/* Add particles randomly to the new vector using the discrete distribution index, higher weighted particles have a stronger
	* likelyhood of being replaced in the new particle vector, this lets the particles that best represent the vehicle location
	* live on. This will also help future particles converge on the correct location and create a smaller distribution that
	* represents the vehicle location most accurately.
	*/
    for (int i = 0; i < particles.size(); i++) {
        resampled_particles.push_back(particles[distribution(gen)]);
		// Reset weight for next updateWeights call
		particles[i].weight = 1.0;
    }

	// Update the particles vector with the resampled particles
    particles = resampled_particles;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
