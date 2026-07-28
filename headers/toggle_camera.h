#ifndef TOGGLE_CAMERA_H
#define TOGGLE_CAMERA_H

#include <glad/glad.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <string>

enum Key_Input {
    W,
    A,
    S,
    D
};

const float YAW         = 180.0f;
const float PITCH       =  0.0f;
const float SPEED       =  2.5f;
const float SENSITIVITY =  0.1f;
const float ZOOM        =  45.0f;

class Camera {
    public:
        glm::vec3 Position;
        glm::vec3 Front;
        glm::vec3 Up;
        glm::vec3 Right;
        glm::vec3 WorldUp;

        float Yaw;
        float Pitch;

        float MovementSpeed;
        float MouseSensitivity;
        float Zoom;

        std::string CurrentView;

        Camera(glm::vec3 position = glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f), float yaw = YAW, float pitch = PITCH) : Front(glm::vec3(0.0f, 0.0f, -1.0f)), MovementSpeed(SPEED), MouseSensitivity(SENSITIVITY), Zoom(ZOOM){
            Position = position;
            WorldUp = up;
            Yaw = yaw;
            Pitch = pitch;
            CurrentView = "front";
            updateCameraVectors();
        }

        void UpdatePosition(glm::vec3 missile_pos){
            if (CurrentView == "front"){
                Position.x = missile_pos.x + 10.0f;
                Position.y = missile_pos.y;
                Position.z = missile_pos.z;
            }
            if (CurrentView == "left"){
                Position.x = missile_pos.x - 3.0f;
                Position.y = missile_pos.y;
                Position.z = missile_pos.z - 12.0f;
            }
            if (CurrentView == "back"){
                Position.x = missile_pos.x - 20.0f;
                Position.y = missile_pos.y;
                Position.z = missile_pos.z;
            }
            if (CurrentView == "right"){
                Position.x = missile_pos.x - 3.0f;
                Position.y = missile_pos.y;
                Position.z = missile_pos.z + 12.0f;
            }
        }

        void ProcessKeyInput(Key_Input input){
            if (input == W){
                CurrentView = "front";
                Yaw = 180.0f;
            }
            if (input == A){
                CurrentView = "left";
                Yaw = 90.0f;
            }
            if (input == S){
                CurrentView = "back";
                Yaw = 0.0f;
            }
            if (input == D){
                CurrentView = "right";
                Yaw = -90.0f;
            }

            Zoom = 45.0f;
            updateCameraVectors();
        }

        glm::mat4 GetViewMatrix(){
            return glm::lookAt(Position, Position + Front, Up);
        }

        void ProcessMouseMovement(float xoffset, float yoffest, GLboolean constrainPitch = true){
            xoffset *= MouseSensitivity;
            yoffest *= MouseSensitivity;

            Yaw += xoffset;
            Pitch += yoffest;

            if (constrainPitch){
                if (Pitch > 89.0f){
                    Pitch = 89.0f;
                }
                if (Pitch < -89.0f){
                    Pitch = -89.0f;
                }
            }

            updateCameraVectors();
        }

        void ProcessMouseScroll(float yoffset){
            Zoom -= 2*(float)yoffset;
            if (Zoom < 1.0f){
                Zoom = 1.0f;
            }
            if (Zoom > 45.0f){
                Zoom = 45.0f;
            }
        }

    private:
        void updateCameraVectors(){
            glm::vec3 front;
            front.x = cos(glm::radians(Yaw)) * cos(glm::radians(Pitch));
            front.y = sin(glm::radians(Pitch));
            front.z = sin(glm::radians(Yaw)) * cos(glm::radians(Pitch));
            Front = glm::normalize(front);

            Right = glm::normalize(glm::cross(Front, WorldUp)); 

            Up    = glm::normalize(glm::cross(Right, Front));
        }
};

#endif