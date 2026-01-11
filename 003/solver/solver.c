#include <stdio.h>
#include <stdlib.h>
// const float one_sixth  = 0x1.555556p-3f; // float 1/6
// const double one_sixth = 0x1.5555555555555p-3; // double 1/6

// include auto-generated lagrangian function:
#include "func_generated.txt"

/*
*
*
*	General solver for N generalized coordinates
*
*
*/

void RK4_N(float* q, float* dq, float* q_dot, float* dq_dot, float t, float dt,
	    void(*dfdt)(float*,float*,float*,float*,float,size_t), size_t N){
	/* RK4 Implementation in N dimensions
	 * q = coordinates array
	 * dq = velocities array
	 * dq_out = derivative of coordinates array
	 * ddq_out = derivative of velocities array
	 * t = current time
	 * dt = time step
	 * dfdt = function that computes derivatives
	 * arguments of dfdt: (q, dq, q_dot, dq_dot, t, N)
	 * N = number of generalized coordinates
	 */

	// Temporary arrays	
	const float one_sixth = 0x1.555556p-3f;
	size_t size = N * sizeof(float);
	float* tmp_q = malloc(size);
	float* tmp_dq = malloc(size);

	// k1, k2, k3, k4 arrays for position and velocity
	float* k1_q = malloc(size); float* k1_dq = malloc(size);
	float* k2_q = malloc(size); float* k2_dq = malloc(size);
	float* k3_q = malloc(size); float* k3_dq = malloc(size);
	float* k4_q = malloc(size); float* k4_dq = malloc(size);

	// Calculate k1, k2, k3, k4
	dfdt(q,dq,k1_q,k1_dq,t,N);
	for(size_t i=0U; i<N; ++i){
		tmp_q[i] = q[i] + 0.5f * dt * k1_q[i];
		tmp_dq[i] = dq[i] + 0.5f * dt * k1_dq[i];
	}
	dfdt(tmp_q,tmp_dq,k2_q,k2_dq,t+0.5f*dt,N);
	for(size_t i=0U; i<N; ++i){
		tmp_q[i] = q[i] + 0.5f * dt * k2_q[i];
		tmp_dq[i] = dq[i] + 0.5f * dt * k2_dq[i];
	}
	dfdt(tmp_q,tmp_dq,k3_q,k3_dq,t+0.5f*dt,N);
	for(size_t i=0U; i<N; ++i){
		tmp_q[i] = q[i] + dt * k3_q[i];
		tmp_dq[i] = dq[i] + dt * k3_dq[i];
	}
	dfdt(tmp_q,tmp_dq,k4_q,k4_dq,t+dt,N);
	// Combine to get final q_dot and dq_dot
	for(size_t i=0U; i<N; ++i){
		q_dot[i] = one_sixth * (k1_q[i] + 2.0f * k2_q[i] + 2.0f * k3_q[i] + k4_q[i]);
		dq_dot[i] = one_sixth * (k1_dq[i] + 2.0f * k2_dq[i] + 2.0f * k3_dq[i] + k4_dq[i]);
	}
	// Cleanup
	free(tmp_q); free(tmp_dq);
	free(k1_q); free(k1_dq);
	free(k2_q); free(k2_dq);
	free(k3_q); free(k3_dq);
	free(k4_q); free(k4_dq);
	return;
}

void next_N(float* q, float* dq, float* new_q, float* new_dq, float dt, size_t N){
	/* Calculating new generalized coordinates */
	int size = N * sizeof(float);
	float* q_dot = malloc(size);
	float* dq_dot = malloc(size);
	RK4_N(q, dq, q_dot, dq_dot, 0.0f, dt, &func, N);
	for (size_t i = 0; i < N; i++)
	{
		new_q[i] = q[i] + q_dot[i] * dt;
		new_dq[i] = dq[i] + dq_dot[i] * dt;
	}
	free(q_dot);
	free(dq_dot);
	return;
}
