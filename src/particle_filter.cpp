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

using std::string;
using std::vector;

static default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles.
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 100;  // TODO: Set the number of particles
  /*
  Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   assuming std[] -> std[std_x, std_y, std_theta]
  */
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  for(int i=0; i<num_particles; ++i){
    Particle p;
    p.id=i;
    //sample from the normal distributions using gen
    p.x=dist_x(gen);
    p.y=dist_y(gen);
    p.theta=dist_theta(gen);
    p.weight=1.0;

    particles.push_back(p);
  }
  is_initialized=true;
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
  for(int j=0; j<num_particles; ++j){
    //Prevent dividing by zero (yaw rate)(Case: when going straight)
    if(abs(yaw_rate<0.001)){
      p[j].x=p[j].x+(velocity*delta_t*cos(p.theta));
      p[j].y=p[j].y+(velocity*delta_t*sin(p.theta));
      p[j].theta= p[j].theta;
    } else{
      //Case: turning
      p[j].x=p[j].x+((velocity/yaw_rate)*(sin(p.[j].theta+(yaw_rate*delta_t))-sin(p.[j].theta)));
      p[j].y=p[j].y+((velocity/yaw_rate)*(cos(p.[j].theta)-cos(p.[j].theta+(yaw_rate*delta_t))));
      p[j].theta=p[j].theta+(yaw_rate*delta_t);
    }
    //Adding noise
    normal_distribution<double> dist_x(p[j].x, std_pos[0]);
    normal_distribution<double> dist_y(p[j].y, std_pos[1]);
    normal_distribution<double> dist_theta(p[j].theta, std_pos[2]);

    p[j].x=dist_x(gen);
    p[j].y=dist_x(gen);
    p[j].x=dist_x(gen);

  }

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
   */
  for(int i=0; i<observations.size(); ++i){
    //Will use the distance formula: distance = sqrt((x2-x1)^2 + (y2=y1)^2)

    //Initial distacne to compare distance value calculated by distance formula and updated.
    double previous_d=numeric<double>std::max();
    double landmark_id=0; //initialize landmark number
    //grab the observation point of interest x and y coordinates
    double observation_x = observations[i].x;
    double observation_y = observations[i].y;

    //Will compare it to every predicted point x and y coordinates
    for(int j=0; j<predicted.size(); ++j){
      double predicted_x = predicted[j].x;
      double predicted_y = predicted[j].y;

      double sub1=observation_x-predicted_x; 
      double sub2=observation_y-predicted_y; 
      double d=sqrt(pow(sub1,2)+pow(sub2,2));

      if(d<previous_d){
        previous_d = d;                 //if we have a new closer prediction, set the distance to newest lowerst distance
        landmark_id = predicted[j].id;  //and set the id to the corresponding predicted id
      }
    }
    observations[i].id = landmark_id;   //set the final lowest distance id to the observed id.
  }

}

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

      Strategy:
        Step1: Convert observations landmarks from car to map coordinates
        Step2: Only use landmarks within sensor range
        Step3: Use kept landmarks in nearest neighboor dataAssociation
        Step4: Calculate the new weights of each partcle
        Step5: Update weights
   */

  double sig_x = std_landmark[0];
  double sig_y = std_landmark[1];

  for(int i=0; i<num_particles; ++i){
    //particles coordinates in map coordinates
    double particle_x=particles[i].x;
    double particle_y=particles[i].y;
    double particle_theta=particles[i].theta;

    vector<LandmarkObs> landmarksInRange; //vector that will contain the landmarks in range

    for(unsigned int k=0; k<map_landmarks.landmarks_list.size(); ++k){
      int landmark_id=map_landmarks[k].id_i;
      double landmark_x=map_landmarks[k].x_f;
      double landmark_y=map_landmarks[k].y_f;

      double d=dist(particle_x, particle_y, landmark_x, landmark_y);

      if(d<sensor_range){
        LandmarkObs landmark;
        landmark.id = landmark_id;
        landmark.x = landmark_x;
        landmark.y = landmark_y;
        landmarksInRange.push_back(landmark);
      }
    }

    vector<LandmarkObs> transformedObs; //vector to contained transformed observations from car to map coordinates

    for(int j=0; j<observations.size(); ++j){
      LandmarkObs transformed_observations;
      transformed_observations.id = j;
      transformed_observations.x = (observations[j].x*cos(particle_theta))-(observations[j].y*sin(particle_theta)) + particle_x;
      transformed_observations.y = (observations[j].x*sin(particle_theta))+(observations[j].y*cos(particle_theta)) + particle_y;
      transformedObs.push_back(transformed_observations);
    }

    //Use nearest neighbor on landmarks in range to associate with observations in the map coordinates
    dataAssociation(landmarksInRange, transformed_observations);

    //Calculate weight of each particle
    particles[i].weight=1.0;

    gauss_norm = 1/(2*M_PI*sig_x*sig_y);

    for(int l=0; l<transformed_observations.size(); ++l){
      double multi_porbability;

      double x_obs = transformed_observations[l].x;
      double y_obs = transformed_observations[l].y;

      for(int m=0; m<landmarksInRange.size(); ++m){
        if(transformed_observations[l].id == landmarksInRange[m].id){
          mu_x = landmarksInRange.x; //coordinates of neerest landmark
          mu_y = landmarksInRange.y; //coordinates of neerest landmark

          double exponent;
          exponent = (pow(x_obs - mu_x, 2) / (2 * pow(sig_x, 2)))
               + (pow(y_obs - mu_y, 2) / (2 * pow(sig_y, 2)));

          multi_porbability = gauss_norm * exp(-exponent);
          particles[i].weight *= multi_porbability;

        }
      }
    }
    weight_normalizer += particles[i].weight
  }
  for(int i=0; i<particles.size(); ++i){
    particles[i].weight/=weight_normalizer;
    weight[i]=particles[i].weight;
  }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

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