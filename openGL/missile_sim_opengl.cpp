#define GLM_ENABLE_EXPERIMENTAL

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>

#include <toggle_camera.h>
#include <shader.h>
#include <normal.h>
#include <csv.h>

class Missile_Obj{
    public:
        unsigned int VBO, VAO;
        glm::vec3 position;
        glm::quat quaternion;
        std::vector<float> delfections;
        float radius;
        float length;

        Missile_Obj(float radius,float length){
            this->VBO = VBO;
            this->VAO = VAO;
            this->radius = radius;
            this->length = length;
            this->position = glm::vec3(0.0f,0.0f,0.0f);
            this->quaternion = glm::quat(1.0f,0.0f,0.0f,0.0f);
            this->delfections = {0.0f,0.0f,0.0f,0.0f};
        }

        std::vector<float> getCircleVerticies(float length){
            std::vector<float> verticies;
            int sections = 30;

            for (float i = 0.0f; i <= sections; i++){
                float theta1 = (i/sections) * (2 * glm::pi<float>());
                float theta2 = ((i+1)/sections) * (2 * glm::pi<float>());

                glm::vec3 v1 = glm::vec3(0.0f,0.0f,-length);
                glm::vec3 v2 = glm::vec3(radius*cos(theta1),radius*sin(theta1),-length);
                glm::vec3 v3 = glm::vec3(radius*cos(theta2),radius*sin(theta2),-length);

                verticies.insert(verticies.end(), {v1.z,v1.y,v1.x});
                verticies.insert(verticies.end(), {v3.z,v3.y,v3.x});
                verticies.insert(verticies.end(), {v2.z,v2.y,v2.x});
            }

            return verticies;
        }

        std::vector<float> getSphereVerticies(){
            std::vector<float> verticies;
            int sections = 30;
            int regions = 30;

            for (float i = 0.0f; i <= sections; i++){
                float theta1 = (i/sections) * glm::pi<float>();
                float theta2 = ((i+1)/sections) * glm::pi<float>();
                for (float j = 0.0f; j < regions; j++){
                    float phi1 = (j/regions) * glm::pi<float>();
                    float phi2 = ((j+1)/regions) * glm::pi<float>();

                    glm::vec3 v1 = glm::vec3(radius*sin(theta1)*cos(phi1),radius*cos(theta1),radius*sin(theta1)*sin(phi1));
                    glm::vec3 v2 = glm::vec3(radius*sin(theta1)*cos(phi2),radius*cos(theta1),radius*sin(theta1)*sin(phi2));
                    glm::vec3 v3 = glm::vec3(radius*sin(theta2)*cos(phi1),radius*cos(theta2),radius*sin(theta2)*sin(phi1));
                    glm::vec3 v4 = glm::vec3(radius*sin(theta2)*cos(phi2),radius*cos(theta2),radius*sin(theta2)*sin(phi2));

                    verticies.insert(verticies.end(), {v1.z,v1.y,v1.x});
                    verticies.insert(verticies.end(), {v2.z,v2.y,v2.x});
                    verticies.insert(verticies.end(), {v3.z,v3.y,v3.x});

                    verticies.insert(verticies.end(), {v2.z,v2.y,v2.x});
                    verticies.insert(verticies.end(), {v4.z,v4.y,v4.x});
                    verticies.insert(verticies.end(), {v3.z,v3.y,v3.x});
                }
            }

            return verticies;
        }

        std::vector<float> getBodyVerticies(){
            std::vector<float> verticies;

            int sections = 30;

            for (float i = 0.0f; i <= sections; i++){
                float theta1 = (i/sections) * (2 * glm::pi<float>());
                float theta2 = ((i+1)/sections) * (2 * glm::pi<float>());

                glm::vec3 v1 = glm::vec3(radius*cos(theta1),radius*sin(theta1),0.0f);
                glm::vec3 v2 = glm::vec3(radius*cos(theta1),radius*sin(theta1),-length);
                glm::vec3 v3 = glm::vec3(radius*cos(theta2),radius*sin(theta2),0.0f);
                glm::vec3 v4 = glm::vec3(radius*cos(theta2),radius*sin(theta2),-length);

                verticies.insert(verticies.end(), {v1.z,v1.y,v1.x});
                verticies.insert(verticies.end(), {v2.z,v2.y,v2.x});
                verticies.insert(verticies.end(), {v3.z,v3.y,v3.x});

                verticies.insert(verticies.end(), {v2.z,v2.y,v2.x});
                verticies.insert(verticies.end(), {v4.z,v4.y,v4.x});
                verticies.insert(verticies.end(), {v3.z,v3.y,v3.x});
            }

            return verticies;
        }

        std::vector<float> getFinVerticies(int num_of_fins){
            std::vector<float> verticies;

            float delta_theta = (2* glm::pi<float>()) / num_of_fins;
            for (int i = 0; i < num_of_fins; i++){
                float theta = 45.0f + delta_theta * i;

                glm::vec3 v1 = glm::vec3(radius*cos(theta+delfections[i]),radius*sin(theta+delfections[i]),-0.5f);
                glm::vec3 v2 = glm::vec3(radius*cos(theta-delfections[i]),radius*sin(theta-delfections[i]),-0.5f-(length/6));
                glm::vec3 v3 = glm::vec3((radius*1.7)*cos(theta-delfections[i]),(radius*1.7)*sin(theta-delfections[i]),-0.5f-(length/6));
                
                glm::vec3 v4 = glm::vec3(radius*cos(theta),radius*sin(theta),-length+(length/6));
                glm::vec3 v5 = glm::vec3(radius*cos(theta),radius*sin(theta),-length);
                glm::vec3 v6 = glm::vec3((radius*1.7)*cos(theta),(radius*1.7)*sin(theta),-length);

                verticies.insert(verticies.end(), {v1.z,v1.y,v1.x});
                verticies.insert(verticies.end(), {v2.z,v2.y,v2.x});
                verticies.insert(verticies.end(), {v3.z,v3.y,v3.x});

                verticies.insert(verticies.end(), {v4.z,v4.y,v4.x});
                verticies.insert(verticies.end(), {v5.z,v5.y,v5.x});
                verticies.insert(verticies.end(), {v6.z,v6.y,v6.x});
            }

            return verticies;
        }

        std::vector<float> getMissileVerticies(){
            std::vector<float> missileVerticies;

            //std::vector<float> frontCapVert = getCircleVerticies(0.0f);
            std::vector<float> backCapVert = getCircleVerticies(length);
            std::vector<float> sphereVert = getSphereVerticies();
            std::vector<float> bodyVert = getBodyVerticies();
            std::vector<float> finVert = getFinVerticies(4);

            missileVerticies.insert(missileVerticies.end(),sphereVert.begin(),sphereVert.end());
            //missileVerticies.insert(missileVerticies.end(),frontCapVert.begin(),frontCapVert.end());
            missileVerticies.insert(missileVerticies.end(),backCapVert.begin(),backCapVert.end());
            missileVerticies.insert(missileVerticies.end(),bodyVert.begin(),bodyVert.end());
            missileVerticies.insert(missileVerticies.end(),finVert.begin(),finVert.end());

            std::vector<float> normals = GetNormals(missileVerticies);
            
            std::vector<float> newVector; 

            for (size_t i = 0; i < missileVerticies.size(); i += 3) {
                newVector.push_back(missileVerticies[i]);
                newVector.push_back(missileVerticies[i + 1]);
                newVector.push_back(missileVerticies[i + 2]);

                newVector.push_back(normals[i]);
                newVector.push_back(normals[i + 1]);
                newVector.push_back(normals[i + 2]);
            }

            missileVerticies = std::move(newVector);
            return missileVerticies;
        }

        void createVertexBufferAndArrays(std::vector<float> verticies){
            glGenBuffers(1, &VBO);
            glGenVertexArrays(1, &VAO);

            glBindBuffer(GL_ARRAY_BUFFER, VBO);
            glBufferData(GL_ARRAY_BUFFER, verticies.size() * sizeof(float), verticies.data(), GL_STATIC_DRAW);

            glBindVertexArray(VAO);

            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
            glEnableVertexAttribArray(0);
            glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
            glEnableVertexAttribArray(1);
        }

        void Draw(std::vector<float> verticies){
            glBindVertexArray(VAO);
            //createVertexBufferAndArrays(verticies);
            glDrawArrays(GL_TRIANGLES, 0, verticies.size());
        }

        std::vector<float> updateParams(int i, std::vector<std::vector<std::string>> data){
            this->position.x = std::stof(data[i][2]);
            this->position.y = std::stof(data[i][4]);
            this->position.z = -std::stof(data[i][3]);

            this->quaternion.w = std::stof(data[i][8]);
            this->quaternion.x = std::stof(data[i][9]);
            this->quaternion.y = std::stof(data[i][11]);
            this->quaternion.z = -std::stof(data[i][10]);

            this->delfections[0] = std::stof(data[i][18]);
            this->delfections[1] = std::stof(data[i][15]);
            this->delfections[2] = std::stof(data[i][16]);
            this->delfections[3] = std::stof(data[i][17]);

            std::vector<float> verticies = getMissileVerticies();
            return verticies;
        }
};

class Target_Obj{
    public:
        unsigned int VAO, VBO;
        float radius;
        glm::vec3 position;

        Target_Obj(float radius){
            this->VBO = VBO;
            this->VAO = VAO;
            this->radius = radius;
            this->position = glm::vec3(0.0f,0.0f,0.0f);
        }

        std::vector<float> GetTargetVerticies(){
            std::vector<float> targetVerticies;

            int sections = 5;
            int regions = 5;

            for (float i = 0.0f; i <= sections; i++){
                float theta1 = (i/sections) * glm::pi<float>();
                float theta2 = ((i+1)/sections) * glm::pi<float>();
                for (float j = 0.0f; j < regions; j++){
                    float phi1 = (j/regions) * (2 * glm::pi<float>());
                    float phi2 = ((j+1)/regions) * (2 * glm::pi<float>());

                    glm::vec3 v1 = glm::vec3(radius*sin(theta1)*cos(phi1),radius*cos(theta1),radius*sin(theta1)*sin(phi1));
                    glm::vec3 v2 = glm::vec3(radius*sin(theta1)*cos(phi2),radius*cos(theta1),radius*sin(theta1)*sin(phi2));
                    glm::vec3 v3 = glm::vec3(radius*sin(theta2)*cos(phi1),radius*cos(theta2),radius*sin(theta2)*sin(phi1));
                    glm::vec3 v4 = glm::vec3(radius*sin(theta2)*cos(phi2),radius*cos(theta2),radius*sin(theta2)*sin(phi2));

                    targetVerticies.insert(targetVerticies.end(), {v1.z,v1.y,v1.x});
                    targetVerticies.insert(targetVerticies.end(), {v2.z,v2.y,v2.x});
                    targetVerticies.insert(targetVerticies.end(), {v3.z,v3.y,v3.x});

                    targetVerticies.insert(targetVerticies.end(), {v2.z,v2.y,v2.x});
                    targetVerticies.insert(targetVerticies.end(), {v4.z,v4.y,v4.x});
                    targetVerticies.insert(targetVerticies.end(), {v3.z,v3.y,v3.x});
                }
            }

            return targetVerticies;
        }

        void createVertexBufferAndArrays(std::vector<float> verticies){
            glGenBuffers(1, &VBO);
            glGenVertexArrays(1, &VAO);

            glBindBuffer(GL_ARRAY_BUFFER, VBO);
            glBufferData(GL_ARRAY_BUFFER, verticies.size() * sizeof(float), verticies.data(), GL_STATIC_DRAW);

            glBindVertexArray(VAO);

            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
            glEnableVertexAttribArray(0);
        }

        void Draw(std::vector<float> verticies){
            glBindVertexArray(VAO);
            glDrawArrays(GL_TRIANGLES, 0, verticies.size());
        }

        void updateParams(int i, std::vector<std::vector<std::string>> data){
            this->position.x = std::stof(data[i][5]);
            this->position.y = std::stof(data[i][7]);
            this->position.z = -std::stof(data[i][6]);
        }
};

unsigned int SCR_WIDTH = 1200;
unsigned int SCR_HEIGHT = 800;

Camera camera(glm::vec3(2.0f, 5000.0f, 10.0f));
float lastX = SCR_WIDTH / 2.0f;
float lastY = SCR_HEIGHT / 2.0f;
bool firstMouse = true;

float deltaTime = 0.0f;	
float lastFrame = 0.0f;

void framebuffer_size_callback(GLFWwindow *window, int width, int height);
void processInput(GLFWwindow *window);
void mouse_callback(GLFWwindow *window, double xpos, double ypos);
void scroll_callback(GLFWwindow *window, double xoffset, double yoffset);

int main(){
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    GLFWwindow *window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "Engagement Simulation", NULL, NULL);
    if (window == NULL){
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)){
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    Shader missileShader("C:\\Software Development\\6 DOF Sim\\shaders\\missile.vs",
                         "C:\\Software Development\\6 DOF Sim\\shaders\\missile.fs");
    Shader targetShader("C:\\Software Development\\6 DOF Sim\\shaders\\target.vs",
                        "C:\\Software Development\\6 DOF Sim\\shaders\\target.fs");

    glViewport(0, 0, SCR_WIDTH, SCR_HEIGHT);

    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
	glfwSetCursorPosCallback(window, mouse_callback);
	glfwSetScrollCallback(window, scroll_callback);

    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

    Missile_Obj missile(0.3f,7.0f);
    std::vector<float> verticies = missile.getMissileVerticies();
    missile.createVertexBufferAndArrays(verticies);
    
    Target_Obj target(5.0f);
    std::vector<float> target_verticies = target.GetTargetVerticies();
    target.createVertexBufferAndArrays(target_verticies);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_FRAMEBUFFER_SRGB);
    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    CSV_Reader reader("C:\\Software Development\\6 DOF Sim\\Data\\GNC_Data2.csv");
    std::vector<std::vector<std::string>> data = reader.getDataVector();
    
    int i = 1;
    double animationTimer = 0;
    double frameDuration = .03;
    while (!glfwWindowShouldClose(window)) {
        float currentFrame = glfwGetTime();
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        animationTimer += deltaTime;

        if (animationTimer >= frameDuration) {
            if (i < data.size() - 1) {
                i += 1;
            } else {
                i = 1;
            }
            animationTimer = 0.0f;
        }

        processInput(window);

        missile.updateParams(i, data);
        target.updateParams(i, data);
        camera.UpdatePosition(missile.position);

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        missileShader.use();

        missileShader.setVec3("viewPos", camera.Position);
        missileShader.setVec3("material.ambient", 0.663f, 0.663f, 0.663f);
		missileShader.setVec3("material.diffuse", 0.663f, 0.663f, 0.663f);
		missileShader.setVec3("material.specular", 0.5f, 0.5f, 0.5f);
		missileShader.setFloat("material.shininess", 32.0f);

        missileShader.setVec3("pointLights[0].position", target.position);
		missileShader.setVec3("pointLights[0].ambient", 0.2f, 0.0f, 0.0f);
		missileShader.setVec3("pointLights[0].diffuse", 0.5f, 0.0f,0.0f);
		missileShader.setVec3("pointLights[0].specular", 1.0f, 1.0f, 1.0f);
		missileShader.setFloat("pointLights[0].constant", 1.0f);
		missileShader.setFloat("pointLights[0].linear", 0.045f);
		missileShader.setFloat("pointLights[0].quadratic", 0.0075f);

        missileShader.setVec3("dirLight.direction", 0.2f, 1.0f, 0.3f);
		missileShader.setVec3("dirLight.ambient", 0.05f, 0.05f, 0.05f);
		missileShader.setVec3("dirLight.diffuse", 0.1f, 0.1f,0.1f);
		missileShader.setVec3("dirLight.specular", 0.1f, 0.1f, 0.1f);

        glm::mat4 projection = glm::perspective(glm::radians(camera.Zoom), (float)SCR_WIDTH/(float)SCR_HEIGHT, 0.1f, 10000.0f);
        glm::mat4 view = camera.GetViewMatrix();
        missileShader.setMat4("projection", projection);
        missileShader.setMat4("view", view);

        glm::mat4 model = glm::translate(glm::mat4(1.0f), missile.position) * glm::toMat4(missile.quaternion);
        missileShader.setMat4("model", model);

        missile.Draw(verticies);

        targetShader.use();
        model = glm::translate(glm::mat4(1.0f), target.position);
        targetShader.setMat4("projection", projection);
        targetShader.setMat4("view", view);
        targetShader.setMat4("model", model);
        target.Draw(target_verticies);

        /* if (i < data.size()-1){
            i += 1;
        }
        else if (i == data.size()-1){
            i = 1;
        } */

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwTerminate();
    return 0;
}

void framebuffer_size_callback(GLFWwindow *window, int width, int height){
	glViewport(0, 0, width, height);
}

void processInput(GLFWwindow *window){
	if(glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
		glfwSetWindowShouldClose(window, true);

	float cameraSpeed = 0.01f + deltaTime;

	if(glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
		camera.ProcessKeyInput(W);

	if(glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
		camera.ProcessKeyInput(S);

	if(glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
		camera.ProcessKeyInput(A);

	if(glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
		camera.ProcessKeyInput(D);
}

void mouse_callback(GLFWwindow *window, double xpos, double ypos){
	if(firstMouse){
		lastX = xpos;
		lastY = ypos;
		firstMouse = false;
	}

	float xoffset = xpos - lastX;
	float yoffest = lastY - ypos;
	lastX = xpos;
	lastY = ypos;

	camera.ProcessMouseMovement(xoffset,yoffest);
}

void scroll_callback(GLFWwindow *window, double xoffset, double yoffest){
	camera.ProcessMouseScroll(yoffest);
}