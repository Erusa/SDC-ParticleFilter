/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"
//#include "debugger_functions.h"
//#include "helper_functions.h"
//#include "multiv_gauss.h"
using std::string;
using std::vector;
using std::normal_distribution;

#define LOG(x) std::cout<<x<<std::endl

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  	std::default_random_engine gen;
  
  	num_particles = 100;  // TODO: Set the number of particles
  	particles.reserve(num_particles);
  	weights.reserve(num_particles);
  	weights.assign(num_particles,1.0);
  	
  // This line creates a normal (Gaussian) distribution
  	normal_distribution<double> dist_x(x, std[0]);
  	normal_distribution<double> dist_y(y, std[1]);
  	normal_distribution<double> dist_theta(theta, std[2]);
  
  // Assign a value to each particle
  	Particle oneParticle;
  	for(int i=0; i<num_particles; ++i){
  		oneParticle.x = dist_x(gen);
    	oneParticle.y = dist_y(gen);
    	oneParticle.theta = dist_theta(gen);
    	oneParticle.weight = 1;
   		oneParticle.id =  i;
    	particles.push_back(oneParticle);
    //std::cout << "Init Particle" <<i<< ":" << particles[i].x << " " << particles[i].y<< " " << particles[i].theta << std::endl;
  		}
	
  	is_initialized= 1;
      // Print your samples to the terminal.
}


void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
    std::default_random_engine gen;
    // This line creates a normal (Gaussian) distribution
  	static normal_distribution<double> dist_x(0, std_pos[0]);
  	static normal_distribution<double> dist_y(0, std_pos[1]);
  	static normal_distribution<double> dist_theta(0, std_pos[2]);
  
  	double x_new, y_new, theta_new;
	for(int i = 0; i<num_particles; ++i){
    	x_new = particles[i].x + (velocity/yaw_rate)*(sin(particles[i].theta+yaw_rate*delta_t) - sin(particles[i].theta)); 
      	x_new += dist_x(gen);
      	particles[i].x =  x_new;
      
    	y_new = particles[i].y + (velocity/yaw_rate)*(-cos(particles[i].theta+yaw_rate*delta_t) + cos(particles[i].theta));
        y_new += dist_y(gen);
      	particles[i].y =  y_new;
      
     	theta_new = particles[i].theta + yaw_rate*delta_t + dist_theta(gen);
      	particles[i].theta = theta_new;
    }
    //std::cout << "Prediction Particle 10 " << particles[10].x << " " << particles[10].y<< " " << particles[10].theta << std::endl;
}

/***************************************************************************************************************/
void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  //COMMENT: Distance from car to landmark is compared. It would be better to compare landmarks positions distance 
  
  vector<LandmarkObs> predictedLandmarks;
  vector<LandmarkObs> observedLandmarks;
  for(int i = 0; i < num_particles; ++i){
    //LOG(i);
    //LOG("predicted Landmarks");
  	predictedLandmarks = getRelevantLandmarks(map_landmarks, sensor_range, particles[i]);
    //LOGLandmarkObsIDs(predictedLandmarks);
    
    //LOG("observed Landmarks");
    observedLandmarks = HomogenousTransformation(observations, particles[i]); //from Car to Map Coordinates
    dataAssociation(predictedLandmarks, observedLandmarks);
    CompleteSetAssociation(particles[i], observedLandmarks);
    //LOGLandmarkObsIDs(observedLandmarks);
    //std::cout<<"association from outside id "<< particles[i].associations[0] <<std::endl;
    //std::cout<<"association from outside x  "<< particles[i].sense_x[0] <<std::endl;
    //std::cout<<"association from outside y "<< particles[i].sense_y[0] <<std::endl;
    //LOG("Updates Weight");
    updateWeight(map_landmarks, observedLandmarks, particles[i], std_landmark);
    //std::cout<<"total weight from hfunction"<< particles[i].weight <<std::endl;
  }
  //LOG("Normalization");
  normalization(); 
  //LOG10ParticlesWeights();
}

vector<LandmarkObs>  ParticleFilter::getRelevantLandmarks(const Map &map_landmarks, const double& sensor_range, const Particle& oneParticle) {
	vector<LandmarkObs> predictedLandmarks;
  	LandmarkObs predictedLandmark;
  	double distance;
   //LOG(map_landmarks.landmark_list.size());
  	for(unsigned int i = 0; i<map_landmarks.landmark_list.size(); ++i ){
    	distance = dist(oneParticle.x, oneParticle.y, map_landmarks.landmark_list[i].x_f, map_landmarks.landmark_list[i].y_f);
      	if (distance < sensor_range){
          predictedLandmark.x = map_landmarks.landmark_list[i].x_f;
          //LOG("predictedLandmark.x");
          //LOG( predictedLandmark.x);
          predictedLandmark.y = map_landmarks.landmark_list[i].y_f;
          predictedLandmark.id =  map_landmarks.landmark_list[i].id_i;  // assumption: id_i = i
          predictedLandmarks.push_back(predictedLandmark);
        }
    }
  return predictedLandmarks;
}      

vector<LandmarkObs> ParticleFilter::HomogenousTransformation(const vector<LandmarkObs> &observationsInCarCoordinates, const Particle& particleInMapCoordinates){
  	
  vector<LandmarkObs> observationsInMapCoordinates;
  observationsInMapCoordinates.reserve(observationsInCarCoordinates.size());
  LandmarkObs observationPoint;
  
  for(unsigned int i = 0; i<observationsInCarCoordinates.size(); ++i){
          observationPoint.x = particleInMapCoordinates.x + cos(particleInMapCoordinates.theta)*observationsInCarCoordinates[i].x - sin(particleInMapCoordinates.theta)*observationsInCarCoordinates[i].y;
          observationPoint.y = particleInMapCoordinates.y + sin(particleInMapCoordinates.theta)*observationsInCarCoordinates[i].x + cos(particleInMapCoordinates.theta)*observationsInCarCoordinates[i].y;                               observationsInMapCoordinates.push_back(observationPoint); 
              //LOG("observationPoint.x");
              //LOG(observationPoint.x);
        }                                                                                 
  
  return  observationsInMapCoordinates;              
  }

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   *  OUTPUT ARGUMENT observation: id was was before zero, id is changed to 
   */
  
  double smallestDistance;
  double Distance;
  double nearstID;
  double observed_x, observed_y;
  for(unsigned int i= 0; i<observations.size(); ++i){
    observed_x = observations[i].x;
    observed_y = observations[i].y;
    
    //init
    smallestDistance = dist(observed_x, observed_y, predicted[0].x, predicted[0].y);
    nearstID = predicted[0].id;
   
  	for(unsigned int j =1; j<predicted.size();++j){
        Distance = dist(observed_x, observed_y, predicted[j].x, predicted[j].y);
    	if(Distance < smallestDistance){
        	nearstID = predicted[j].id;
          	smallestDistance = Distance;
        }    
    }
    observations[i].id = nearstID;
    //std::cout<<"obs_id from" <<i<<": "<< observations[i].id;
  } 
  //std::cout<<std::endl;
}

void ParticleFilter::updateWeight(const Map &map_landmarks, const vector<LandmarkObs>& observedLandmarks, Particle &oneParticle, const double std_landmark[]){
  oneParticle.weight =1.0;
  int landMarkID;
  double predicted_x;
  double predicted_y;
  for(unsigned int i = 0; i<observedLandmarks.size(); ++i){
     landMarkID = observedLandmarks[i].id-1;
     predicted_x = map_landmarks.landmark_list[landMarkID].x_f; 
     predicted_y = map_landmarks.landmark_list[landMarkID].y_f; 
  	 oneParticle.weight *= multiv_prob(std_landmark[0], std_landmark[1], observedLandmarks[i].x, observedLandmarks[i].y, predicted_x, predicted_y);
     //std::cout<< "weight "<< i<<": "<< oneParticle.weight<<"--> p_x: "<<predicted_x << " p_y: "<<predicted_y << " o_x: "<< observedLandmarks[i].x << " o_y: " <<observedLandmarks[i].y << std::endl;
  }
     //std::cout<< "total weight "<< ": "<< oneParticle.weight<<std::endl;
}

void ParticleFilter::normalization(){
	double sum = 0;
  for(int i = 0; i<num_particles; ++i){
    	sum +=  particles[i].weight;
    }
  for(int i = 0; i<num_particles; ++i){
    	particles[i].weight =  particles[i].weight / sum;
    }
}

/***************************************************************************************************************/
void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  std::default_random_engine generator;
  //LOG("Resample");
  vector<int> particlesProbabilities;
  particlesProbabilities.reserve(num_particles);
  setParticlesProbabilities(particlesProbabilities); // here discrete distribution is set using particles[i].weight TODO!!!!!
  std::discrete_distribution<int> distribution(particlesProbabilities.begin(), particlesProbabilities.end());
  
  for(int i =0; i<num_particles; ++i){
    int random_ID = distribution(generator); // here one particule is randomly choseen
     //std::cout<< "old weight "<< i<<": "<< weights[i]<<std::endl;
  	weights[i] = random_ID;
    //std::cout<< "weight "<< i<<": "<< weights[i]<<std::endl;
  }

  actualizeParticles(); // aqui se rescriba nuevas particulas donde varias se duplican
}

void ParticleFilter::actualizeParticles(){
    // Assign a value to each particle
  Particle oneParticle;
   vector<Particle> newParticles;
   newParticles.reserve(num_particles);
  int ParticleID;
  for(int i=0; i<num_particles; ++i){
    ParticleID = weights[i];
  	oneParticle.x = particles[ParticleID].x;
    oneParticle.y = particles[ParticleID].y;
    oneParticle.theta = particles[ParticleID].theta;
    oneParticle.weight = particles[ParticleID].weight;
    oneParticle.id =  i;
    newParticles.push_back(oneParticle);
    //std::cout << "Init Particle" <<i<< ":" << particles[i].x << " " << particles[i].y<< " " << particles[i].theta << std::endl;
  }
	
  particles = newParticles;
}

void ParticleFilter::setParticlesProbabilities(vector<int> &particlesProbabilities){
  for(int i=0; i<num_particles; ++i){
  	particlesProbabilities.push_back((int)(round(particles[i].weight * 1000000)));}
}
  
void ParticleFilter::LOG10ParticlesWeights(){
	std::cout<< "Printing Particle Weights" <<std::endl;
  for(int i = 0; i<num_particles; i=i+10){
    std::cout << "Particle"<< i << ":" << particles[i].weight <<std::endl;
  }
}
/***************************************************************************************************************/
void ParticleFilter::CompleteSetAssociation(Particle& particle, const vector<LandmarkObs> &observedLandmarks)
{
	vector<int> associations;
    vector<double> sense_x;
    vector<double> sense_y;
    unsigned int n = observedLandmarks.size();
  	associations.reserve(n);
  	sense_x.reserve(n);
  	sense_y.reserve(n);
  
  for(unsigned int i=0; i<n; ++i){
  	sense_x.push_back(observedLandmarks[i].x);
    sense_y.push_back(observedLandmarks[i].y);
    associations.push_back(observedLandmarks[i].id);
    //std::cout<<"asc_id" <<i<<": "<< associations[i];
  }
  
  SetAssociations(particle, associations, sense_x, sense_y);
   //std::cout<<std::endl;
}


void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}
string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}