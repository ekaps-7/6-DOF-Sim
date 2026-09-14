#include <IGNC.cpp>

int main(int argc,char* argv[]){
    if (argc < 2)
    {
        std::cerr << "Enter a cfg file\n";
        return 1;
    }

    Simulation sim;
    
    sim.run(argv[1]);

/* 
    sim.run_Monte_Carlo(StateModel::LR_HS_LATERAL_CROSSING,
                        StateModel::FULL_BODY,
                        100);
 */
    return 0;
}