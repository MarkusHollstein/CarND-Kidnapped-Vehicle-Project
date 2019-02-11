/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
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
  	
	num_particles = 50;//don't know yet
  	//create normal distributions for x,y,theta
  	default_random_engine gen;
  	normal_distribution<double> dist_x(x, std[0]);
  	normal_distribution<double> dist_y(y, std[1]);
  	normal_distribution<double> dist_theta(theta, std[2]);
  	for(unsigned int i= 0;i<num_particles;i++)
    {
      double part_x = dist_x(gen);
      double part_y = dist_y(gen);
      
      double part_theta = fmod(dist_theta(gen),2*M_PI);
      //initialize new particle and assign values
      Particle part;
      part.x = part_x;
      part.y = part_y;
      part.theta = part_theta;
      part.weight=1;
      part.id = i;
      particles.push_back(part);
      weights.push_back(1.);
    }
  	is_initialized =true; 
  	
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
  	

  	for(unsigned int i=0;i<num_particles;i++)
    {
      	Particle part = particles[i];  
      	if(yaw_rate==0)
        {
          	part.x +=delta_t*velocity*cos(part.theta);          	
          	part.y +=delta_t*velocity*sin(part.theta);
          	//theta remains unchanged
        }
      	else
        {      		
    		part.x += velocity/yaw_rate *(sin(part.theta+yaw_rate*delta_t)-sin(part.theta));
       		part.y += velocity/yaw_rate * (cos(part.theta)-cos(part.theta+yaw_rate*delta_t));
      		part.theta += yaw_rate*delta_t;   
        }
      	std::random_device rd{};
    	std::mt19937 gen{rd()};
  		normal_distribution<double> dist_x{part.x, std_pos[0]};
  		normal_distribution<double> dist_y(part.y, std_pos[1]);
  		normal_distribution<double> dist_theta(part.theta, std_pos[2]);
      	part.x = dist_x(gen); 
      	part.y = dist_y(gen);
      	part.theta = dist_theta(gen);
      	//part.theta = fmod(part.theta,2*M_PI);
      	particles[i] = part;
    }
	
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
  	
  	for(unsigned int i= 0;i<observations.size();i++)//
    {       
      double min_dist=dist(observations[i].x,observations[i].y,predicted[0].x,predicted[0].y);
      int min_idx=0;      
      for(unsigned int j=1;j<predicted.size();j++)
      {
        double curr_dist = dist(observations[i].x,observations[i].y,predicted[j].x,predicted[j].y);
        if(curr_dist < min_dist )
        {
          min_dist = curr_dist;
          min_idx = j;
        }        
      }
      observations[i].id = min_idx;     
    }	
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
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
  
  //plan: coordinate transformation; data association function(transf_landmarks,observations) to get the associated landmarks
  // to the observations; finally calculate exp stuff...  
  for(unsigned int k=0;k<num_particles;k++)
  {
  	Particle part=particles[k]; 
  	//transform observations to (global) map coordinates
  	std::vector<LandmarkObs> transformed;
  	for(unsigned int i=0;i<observations.size();i++)
  	{
      	double x_map = cos(part.theta)*observations[i].x -sin(part.theta)*observations[i].y + part.x;
    	double y_map = sin(part.theta)*observations[i].x +cos(part.theta)*observations[i].y +part.y;      	
   		LandmarkObs curr;
    	curr.x = x_map;
    	curr.y = y_map;
      	curr.id = -1;
    	transformed.push_back(curr);
    	//id will be determined by association
  	}  	  
  	std::vector<LandmarkObs> landmarks;    
  	for(unsigned int i=0;i<map_landmarks.landmark_list.size();i++)
  	{
      
    	LandmarkObs curr;
    	curr.x = map_landmarks.landmark_list[i].x_f;
    	curr.y = map_landmarks.landmark_list[i].y_f;
   	 	curr.id = map_landmarks.landmark_list[i].id_i;  
      	if(dist(curr.x,curr.y,part.x,part.y)<=1*sensor_range)
        {
    		landmarks.push_back(curr);
        }      	
  	}
  	dataAssociation(landmarks,transformed);
  	double weight = 1;
  	for(unsigned int i= 0 ; i< transformed.size();i++)
  	{      
      double expo = -(pow(landmarks[transformed[i].id].x-transformed[i].x,2)/(pow(std_landmark[0],2))+pow(landmarks[transformed[i].id].y-transformed[i].y,2)/(pow(std_landmark[0],2)))/2;      
      weight *= exp(expo)/(2*M_PI*std_landmark[0]*std_landmark[1]);  	  
    }  
    particles[k].weight = weight;
    weights[k]=weight;
  }  
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
  	auto it_begin =weights.begin();
  	auto it_end = weights.end();  	
    std::random_device rd;
    std::mt19937 gen(rd());
  	std::discrete_distribution<> dist(it_begin,it_end);
  	std::vector<Particle> new_particles;      
	for(unsigned int i=0; i<num_particles;i++)
    {
      int idx = dist(gen);
      Particle part=particles[idx];
      part.id = i;
      new_particles.push_back(part);
      weights[i]=part.weight;
    }
  	particles = new_particles;  	
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
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
