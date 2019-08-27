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


using namespace std;  

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

  //Creates normal Gaussian distribution for x, y, and theta
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  //For every particle filter, assign an id, x, y, theta, and weight and add them
  //to the particles 
  for(int i=0; i<num_particles; ++i){
    Particle p;                //Initialize particle p to type Particle
    p.id=i;                    //set the id of the particle to the particle number as it was created starting at id=0
    //std::cout<<"particle id: "<<i<<std::endl;
    //sample from the normal distributions using gen using the provided gps data (x, y, theta, std)
    p.x=dist_x(gen);           //set x coordinate for this particle randomly in (map coordinates? Probably)
    p.y=dist_y(gen);           //set y coordinate for this particle randomly in (map coordinates? Probably)
    p.theta=dist_theta(gen);   //set theta for this particle randomly in (map coordinates? Probably)
    p.weight=1.0;              //set every particle weight to 1.0

    particles.push_back(p);    //add the partcile to the Particle type vector "particles"
  }
  //std::cout<<"Number of particles: "<<particles.size()<<std::endl;
  is_initialized=true;         //set initialization to done (true)

  //std::cout<<"Initialization complete!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
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

  //For every particle, update x, y, and yaw angle using the prediciton bicycle model
  //std::cout<<"Starting prediction"<<std::endl;
  for(int j=0; j<num_particles; ++j){
    //std::cout<<"Particle: "<<j<<std::endl;
    //Prevent dividing by zero (yaw rate)(Case: when going straight)
    if(abs(yaw_rate)<0.00001){
      //std::cout<<"Going straight"<<std::endl;
      particles[j].x=particles[j].x+(velocity*delta_t*cos(particles[j].theta));
      particles[j].y=particles[j].y+(velocity*delta_t*sin(particles[j].theta));
      particles[j].theta= particles[j].theta;
    } else{
      //Case: turning
      //std::cout<<"Turning"<<std::endl;
      particles[j].x=particles[j].x+((velocity/yaw_rate)*(sin(particles[j].theta+(yaw_rate*delta_t))-sin(particles[j].theta)));
      particles[j].y=particles[j].y+((velocity/yaw_rate)*(cos(particles[j].theta)-cos(particles[j].theta+(yaw_rate*delta_t))));
      particles[j].theta=particles[j].theta+(yaw_rate*delta_t);
    }
    //Adding noise
    normal_distribution<double> dist_x(particles[j].x, std_pos[0]);
    normal_distribution<double> dist_y(particles[j].y, std_pos[1]);
    normal_distribution<double> dist_theta(particles[j].theta, std_pos[2]);

    particles[j].x=dist_x(gen);
    particles[j].y=dist_y(gen);
    particles[j].theta=dist_theta(gen);
  }
  //std::cout<<"Prediction complete!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
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
  //std::cout<<"STARTING DATA ASSOCIATION!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;

  for(int i=0; i<observations.size(); ++i){
    //Will use the distance formula: distance = sqrt((x2-x1)^2 + (y2=y1)^2)

    //Initial distacne to compare distance value calculated by distance formula and updated.
    double previous_d=numeric_limits<double>::max();              //Initialize previous d to super high number initially
    double landmark_id=0;                                         //initialize landmark number
    //grab the observation point of interest x and y coordinates (Should be in map coordinates now)
    double observation_x = observations[i].x;
    double observation_y = observations[i].y;

    //Will compare it to every predicted point x and y coordinates
    for(int j=0; j<predicted.size(); ++j){
      double predicted_x = predicted[j].x;
      double predicted_y = predicted[j].y;

      double d=dist(observation_x, observation_y, predicted_x, predicted_y);

      if(d<previous_d){
        previous_d = d;                 //if we have a new closer prediction, set the distance to newest lowerst distance
        landmark_id = predicted[j].id;  //and set the id to the corresponding predicted id
      }
    }
    observations[i].id = landmark_id;   //set the final lowest distance id to the observed id.
    //std::cout<<"New Observation IDs: "<<observations[i].id<<std::endl;
  }
  //std::cout<<"ENDING DATA ASSOCIATION!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
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
        
  double sig_x = std_landmark[0];           //set the std of x to sig_x
  //std::cout<<"sig_x: "<<sig_x<<std::endl;
  double sig_y = std_landmark[1];           //set the std of y to sig_y
  //std::cout<<"sig_y: "<<sig_y<<std::endl;

  double gauss_norm = 1/(2*M_PI*sig_x*sig_y);
  //gauss_norm constant 1.76839

  //For every particle
  for(int i=0; i<num_particles; ++i){
    //std::cout<<"particle id: "<<i<<std::endl;
    //particles coordinates in map coordinates
    double multi_porbability;

    //Grab the particles x, y, and theta in vehicle coordinates
    double particle_x=particles[i].x;
    double particle_y=particles[i].y;
    double particle_theta=particles[i].theta;

    //std::cout<<"------particle (x, y, theta): "<<particle_x<<" "<<particle_y<<" "<<particle_theta<<std::endl;

    //-------------------------------------Work only with the landmarks in range of hte sensor----------------------------------------
    vector<LandmarkObs> landmarksInRange; //vector that will contain the landmarks in range

    for(unsigned int k=0; k<map_landmarks.landmark_list.size(); ++k){
      //std::cout<<"------map landmark "<<k<<std::endl;

      int landmark_id=map_landmarks.landmark_list[k].id_i;
      double landmark_x=map_landmarks.landmark_list[k].x_f;
      double landmark_y=map_landmarks.landmark_list[k].y_f;

      //std::cout<<"------landmark (landmark id, landmark_x, landmark_y): "<<landmark_id<<" "<<landmark_x<<" "<<landmark_y<<std::endl;

      double d=dist(particle_x, particle_y, landmark_x, landmark_y);
      
      //std::cout<<"------distance calculated between particle and lanmark 'd' vs. sensor range: "<<d<<" vs. "<<sensor_range<<std::endl;      

      if(d<=sensor_range){
        LandmarkObs landmark;
        landmark.id = landmark_id;
        //std::cout<<"Landmark within sensor range: "<<landmark.id<<std::endl;
        landmark.x = landmark_x;
        landmark.y = landmark_y;
        landmarksInRange.push_back(landmark);
        //std::cout<<"Length of landmarks in range "<<landmarksInRange.size()<<std::endl;
      }
    }
    //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    //------------------------------Change landmark observations from car coordinates to map coordinates---------------------
    vector<LandmarkObs> transformedObs; //vector to contained transformed observations from car to map coordinates

    for(int j=0; j<observations.size(); ++j){
      LandmarkObs transformed_observations;
      transformed_observations.id = j;
      transformed_observations.x = (observations[j].x*cos(particle_theta))-(observations[j].y*sin(particle_theta)) + particle_x;
      transformed_observations.y = (observations[j].x*sin(particle_theta))+(observations[j].y*cos(particle_theta)) + particle_y;
      transformedObs.push_back(transformed_observations);
    }
    //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    //---------------------Use nearest neighbor on landmarks in range to associate with observations in the map coordinates-------------
    dataAssociation(landmarksInRange, transformedObs);
    // transformedObs should now contain the observations in map coordinates with thet id of their cooresponding nearest landmark
    //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    //-----------------------------------------------Calculate weight of each particle----------------------------------------------------
    particles[i].weight=1.0; //reset the weight of each particle

    //Iterate through each observation
    for(int l=0; l<transformedObs.size(); ++l){ 
      //grab current obesrvation coordinates
      double x_obs = transformedObs[l].x;
      double y_obs = transformedObs[l].y;
      
      //Iterate through each landmark in range
      for(int m=0; m<landmarksInRange.size(); ++m){
        //Once you find the matching ID's compute the weight or "multi_porbability"
        if(transformedObs[l].id == landmarksInRange[m].id){
          double mu_x;
          double mu_y;

          mu_x = landmarksInRange[m].x; //coordinates of nearest landmark
          mu_y = landmarksInRange[m].y; //coordinates of nearest landmark

          double inExponent = ((pow(x_obs - mu_x, 2) / (2 * pow(sig_x, 2)))
               + (pow(y_obs - mu_y, 2) / (2 * pow(sig_y, 2))));
          std::cout<<"inExponent: "<<inExponent<<std::endl;
          std::cout<<"Exponenet: "<<exp(-inExponent)<<std::endl;
          multi_porbability = gauss_norm * exp(-inExponent);
        }//End of finding mathc
      }//End of going through each lanmark in range

      //std::cout<<"--------------weight: "<<multi_porbability<<std::endl;
      particles[i].weight *= multi_porbability;
    }//End of going through each observation
    //std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Particle "<<i<<" weight "<<particles[i].weight<<std::endl;
  }//End of each particle
  //std::cout<<"ENED OF UPDATE WEIGHTS!!!"<<std::endl;
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  
  vector<Particle> resampled_particles;

  vector<double> weights;
  for(int i =0; i<num_particles; ++i){
    weights.push_back(particles[i].weight);
  }

  //random particle index
  uniform_int_distribution<int> particle_index(0, num_particles-1);
  default_random_engine gen;

  double beta = 0.0;

  double max_weight = 2.0 * *max_element(weights.begin(), weights.end());

  int index = particle_index(gen);

  uniform_real_distribution<double> random_weight(0.0, max_weight);

  for(int i=0; i<num_particles; ++i){
    beta += random_weight(gen) * max_weight;

    while(beta> weights[index]){
      beta -= weights[index];
      index=(index+1) % num_particles;
    }
    resampled_particles.push_back(particles[index]);
  }  
  particles = resampled_particles;
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

  //particle.associations.clear();
  //particle.sense_x.clear();
  //particle.sense_y.clear();
  

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