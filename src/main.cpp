#include <uWS/uWS.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "helpers.h"
#include "json.hpp"
#include "spline.h"

// for convenience
using nlohmann::json;
using std::string;
using std::vector;

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  std::ifstream in_map_(map_file_.c_str(), std::ifstream::in);

  string line;
  while (getline(in_map_, line)) {
    std::istringstream iss(line);
    double x;
    double y;
    float s;
    float d_x;
    float d_y;
    iss >> x;
    iss >> y;
    iss >> s;
    iss >> d_x;
    iss >> d_y;
    map_waypoints_x.push_back(x);
    map_waypoints_y.push_back(y);
    map_waypoints_s.push_back(s);
    map_waypoints_dx.push_back(d_x);
    map_waypoints_dy.push_back(d_y);
  }
  
  int lane = 1; //initial lane
  int last_lane = 2;
  double ref_vel = 0.0; // miles per hour
  
  string state = "KL";
  
  // max_vel = 50 miles/h
  // max_acc = 10 m/s2
  // max_jerk = 10 m/s3
  
  h.onMessage([&ref_vel,&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,
               &map_waypoints_dx,&map_waypoints_dy,&lane,&state,&last_lane]
              (uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
               uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
          // Main car's localization Data
          double car_x = j[1]["x"];
          double car_y = j[1]["y"];
          double car_s = j[1]["s"];
          double car_d = j[1]["d"];
          double car_yaw = j[1]["yaw"];
          double car_speed = j[1]["speed"];

          // Previous path data given to the Planner
          auto previous_path_x = j[1]["previous_path_x"];
          auto previous_path_y = j[1]["previous_path_y"];
          // Previous path's end s and d values 
          double end_path_s = j[1]["end_path_s"];
          double end_path_d = j[1]["end_path_d"];

          // Sensor Fusion Data, a list of all other cars on the same side 
          //   of the road.
          auto sensor_fusion = j[1]["sensor_fusion"];

          json msgJson;

          /**
           * TODO: define a path made up of (x,y) points that the car will visit
           *   sequentially every .02 seconds
           */
                           
          //define the actual (x,y) points we will use for the planner
          vector<double> next_x_vals;
          vector<double> next_y_vals;

          int prev_size = previous_path_x.size();
          
          if (prev_size>0){
            car_s = end_path_s;
          }
          
          bool too_close = false;
          double vx;
          double vy;
          double check_speed;
          double check_car_s;
          float d;
          //find ref_v to use
          for (int i=0;i<sensor_fusion.size();i++){
            //check if another car is in my car's lane
            d = sensor_fusion[i][6];
            if ((d<(2+4*lane+2)) && (d>(2+4*lane-2))){
              vx = sensor_fusion[i][3];
              vy = sensor_fusion[i][4];
              check_speed = sqrt(vx*vx + vy*vy);
              check_car_s = sensor_fusion[i][5];
              
              //if using previous points we can project s value outwards in time, meaning that we can predict where the obstacle car will be next 
              check_car_s += (double)prev_size*.02*check_speed;
              //check if the other car is in front of us and whether it is too close e.g. 30m
              if ((check_car_s > car_s) && ((check_car_s - car_s)<30)){
                too_close = true;
                
                //basic logic to change to a faster lane. TODO: use path planner instead
                if (lane>0){
                  lane = 0;
                }
              }
            }
          }
          
          
          double COLLISION_FACTOR = 3.0;
          double SPEED_FACTOR = 1.0;
          double SAFETY_DISTANCE = 25.0;
          double LOOK_AHEAD_TIME = 0.04;
          
          //simple cost function planner
          double next_lane_cost;
          double previous_lane_cost;
          double KL_cost = SPEED_FACTOR*speed_cost(ref_vel);
          if (too_close){
            vector<vector<double>> predictions = get_predictions(prev_size,sensor_fusion,LOOK_AHEAD_TIME);
            ref_vel -= 0.224; //approx -5 meters/s^2
            
            //check if it's possible to change to the lane in the right
            if (lane == 0){
              next_lane_cost = check_next_lane(lane,car_s,predictions,COLLISION_FACTOR,SAFETY_DISTANCE);
              if (next_lane_cost < KL_cost){
                lane++;
              }
            }

            // check if either adjacent lane is better than the current one
            else if ((lane>0)&&(lane < last_lane)){
              previous_lane_cost = check_prev_lane(lane,car_s,predictions,COLLISION_FACTOR,SAFETY_DISTANCE);
              next_lane_cost = check_next_lane(lane,car_s,predictions,COLLISION_FACTOR,SAFETY_DISTANCE);

              if ( (previous_lane_cost < next_lane_cost) && (previous_lane_cost < KL_cost) ){
                lane--;
              }
              else if (next_lane_cost < KL_cost){
                lane++;
              }
            }
            else //lane == last_lane. //check if it's possible to change to the lane in the left
            {
              previous_lane_cost = check_prev_lane(lane,car_s,predictions,COLLISION_FACTOR,SAFETY_DISTANCE);
              if (previous_lane_cost < KL_cost){
                lane--;
              }
            }
            
          }
          else if (ref_vel < 49.5){
            ref_vel += 0.224; //approx +5 meters/s^2
          }  
                  
          /*
          if (too_close){
            ref_vel -= 0.224; //approx -5 meters/s^2
          }
          else if (ref_vel < 49.5){
            ref_vel += 0.224; //approx +5 meters/s^2
          }*/
          
                  
          // create a list of widely spaced (x,y) waypoints, evenly spaced at 30m
          vector<double> pts_x;
          vector<double> pts_y;
          
          //reference x, y, and yaw states. Either we will reference the starting point as where the car is or at the previous path's end point
          double ref_x = car_x;
          double ref_y = car_y;
          double ref_yaw = deg2rad(car_yaw);
          
          // use 2 points that make the path tangent to the previous path's end point
          if (prev_size < 2){
            double prev_car_x = car_x - cos(car_yaw);
            double prev_car_y = car_y - sin(car_yaw);
            pts_x.push_back(prev_car_x);
            pts_x.push_back(car_x);
            pts_y.push_back(prev_car_y);
            pts_y.push_back(car_y);
          }
          else
          {
            ref_x = previous_path_x[prev_size - 1];
            ref_y = previous_path_y[prev_size - 1];
            
            double ref_x_prev = previous_path_x[prev_size - 2];
            double ref_y_prev = previous_path_y[prev_size - 2];
            ref_yaw = atan2(ref_y - ref_y_prev,ref_x - ref_x_prev);
            
            pts_x.push_back(ref_x_prev);
            pts_x.push_back(ref_x);
            pts_y.push_back(ref_y_prev);
            pts_y.push_back(ref_y);
          }
          
          // In Frenet add evenly 30m spaced points ahead of the starting reference
          vector<double> next_wp0 = getXY(car_s+30,(2+4*lane),map_waypoints_s,map_waypoints_x,map_waypoints_y);
          vector<double> next_wp1 = getXY(car_s+60,(2+4*lane),map_waypoints_s,map_waypoints_x,map_waypoints_y);
          vector<double> next_wp2 = getXY(car_s+90,(2+4*lane),map_waypoints_s,map_waypoints_x,map_waypoints_y);
          
          pts_x.push_back(next_wp0[0]);
          pts_x.push_back(next_wp1[0]);
          pts_x.push_back(next_wp2[0]);
          
          pts_y.push_back(next_wp0[1]);
          pts_y.push_back(next_wp1[1]);
          pts_y.push_back(next_wp2[1]);
          
          //shift car reference angle to 0 degrees
          for (int i=0;i<pts_x.size();i++){
            double shift_x = pts_x[i]-ref_x;
            double shift_y = pts_y[i]-ref_y;
            
            pts_x[i] = (shift_x*cos(0-ref_yaw)-shift_y*sin(0-ref_yaw));
            pts_y[i] = (shift_x*sin(0-ref_yaw)+shift_y*cos(0-ref_yaw));
          }
          
          //create a splne
          tk::spline s;
          
          // set (x,y) points to the spline
          s.set_points(pts_x,pts_y);
          
          //start with all of the previous path points from last time
          for (int i = 0; i < previous_path_x.size(); i++){
            next_x_vals.push_back(previous_path_x[i]);
            next_y_vals.push_back(previous_path_y[i]);
          }
          
          //calculate how to break up spline points so that we can travel at our desired reference velocity
          double target_x = 30.0;
          double target_y = s(target_x);
          double target_dist = sqrt((target_x)*(target_x) + (target_y)*(target_y));
          
          double x_add_on = 0.0;
          
          //fill up the rest of our path planner after filling it with previous points. Here we will always output 50 points
          for (int i = 1; i <= 50-previous_path_x.size(); i++){
            double N = (target_dist/(.02*ref_vel/2.24));
            double x_point = x_add_on + (target_x)/N;
            double y_point = s(x_point);
            
            x_add_on = x_point;
            
            double x_ref = x_point;
            double y_ref = y_point;
            
            //rotate back to world coordinates
            x_point = (x_ref*cos(ref_yaw)-y_ref*sin(ref_yaw));
            y_point = (x_ref*sin(ref_yaw)+y_ref*cos(ref_yaw));
            
            x_point += ref_x;
            y_point += ref_y;
            
            next_x_vals.push_back(x_point);
            next_y_vals.push_back(y_point);
            
          }
          
          /*
          // Stay in the center of the lane
          double dist_inc = 0.5;
          for (int i = 0; i < 50; ++i) {
            double next_s = car_s+((i+1)*dist_inc);
            double next_d = 6; //stay in the center of the middle lane
            vector<double> xy = getXY(next_s,next_d,map_waypoints_s,map_waypoints_x,map_waypoints_y);
            next_x_vals.push_back(xy[0]);
            next_y_vals.push_back(xy[1]);
          }
          */

          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;

          auto msg = "42[\"control\","+ msgJson.dump()+"]";

          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }  // end "telemetry" if
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }  // end websocket if
  }); // end h.onMessage

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  
  h.run();
}