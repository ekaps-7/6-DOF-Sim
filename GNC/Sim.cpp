#include <IGNC.cpp>

int main(){
    Simulation sim;
    
    sim.run(StateModel::LR_HS_LATERAL_CROSSING,
            StateModel::TRIM);

/* 
    sim.run_Monte_Carlo(StateModel::LR_HS_LATERAL_CROSSING,
                        StateModel::FULL_BODY,
                        100);
 */
    return 0;
}