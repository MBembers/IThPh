#include <stdio.h>
#include <stdlib.h>
// const float one_sixth  = 0x1.555556p-3f; // float 1/6
// const double one_sixth = 0x1.5555555555555p-3; // double 1/6

// dfdt func:
/*
 * Auto-generated function to compute derivatives
 */

// Simple pendulum
void func(float* q, float* dq, float* _dq, float* _ddq, float t, size_t N) {
    // Auto-generated Euler-Lagrange Equations using sympy.physics.mechanics
    // Constants have been collapsed into their values.
    _dq[0] = dq[0];
    _ddq[0] = -9.8100000000000005*sin(q[0]);
return;
}


/* --- 1D Functions ---
                                                                                          
   ▄▄▄     ▄▄▄▄▄                                                                          
  █▀██     ██▀▀▀██                                             ██                         
    ██     ██    ██            ▄▄█████▄  ▀██  ███  ▄▄█████▄  ███████    ▄████▄   ████▄██▄ 
    ██     ██    ██            ██▄▄▄▄ ▀   ██▄ ██   ██▄▄▄▄ ▀    ██      ██▄▄▄▄██  ██ ██ ██ 
    ██     ██    ██             ▀▀▀▀██▄    ████▀    ▀▀▀▀██▄    ██      ██▀▀▀▀▀▀  ██ ██ ██ 
 ▄▄▄██▄▄▄  ██▄▄▄██             █▄▄▄▄▄██     ███    █▄▄▄▄▄██    ██▄▄▄   ▀██▄▄▄▄█  ██ ██ ██ 
 ▀▀▀▀▀▀▀▀  ▀▀▀▀▀                ▀▀▀▀▀▀      ██      ▀▀▀▀▀▀      ▀▀▀▀     ▀▀▀▀▀   ▀▀ ▀▀ ▀▀ 
                                          ███                                             
                                                                                          
 */
void RK4_1D(float* x, float* v, float* dx, float* dv, float t, float dt,
	    void(*dfdx)(float*,float*,float*,float*,float,size_t), size_t N){
	/* RK4 Implementation in 1D
	 * x = position array
	 * v = velocity array
	 * dx = derivative of position array
	 * dv = derivative of velocity array
	 * t = current time
	 * dt = time step
	 * dfdx = function that computes derivatives
	 * arguments of dfdx: (x, v, dx, dv, t, N)
	 * N = number of elements
	 */

	// Temporary arrays	
	const float one_sixth = 0x1.555556p-3f;
	size_t size = N * sizeof(float);
	float* tmp_x = malloc(size);
	float* tmp_v = malloc(size);

	// k1, k2, k3, k4 arrays for position and velocity
	float* k1_dx = malloc(size); float* k1_dv = malloc(size);
	float* k2_dx = malloc(size); float* k2_dv = malloc(size);
	float* k3_dx = malloc(size); float* k3_dv = malloc(size);
	float* k4_dx = malloc(size); float* k4_dv = malloc(size);

	// Calculate k1, k2, k3, k4
	dfdx(x,v,k1_dx,k1_dv,t,N);
	for(size_t i=0U; i<N; ++i){
		tmp_x[i] = x[i] + 0.5f * dt * k1_dx[i];
		tmp_v[i] = v[i] + 0.5f * dt * k1_dv[i];
	}
	dfdx(tmp_x,tmp_v,k2_dx,k2_dv,t+0.5f*dt,N);
	for(size_t i=0U; i<N; ++i){
		tmp_x[i] = x[i] + 0.5f * dt * k2_dx[i];
		tmp_v[i] = v[i] + 0.5f * dt * k2_dv[i];
	}
	dfdx(tmp_x,tmp_v,k3_dx,k3_dv,t+0.5f*dt,N);
	for(size_t i=0U; i<N; ++i){
		tmp_x[i] = x[i] + dt * k3_dx[i];
		tmp_v[i] = v[i] + dt * k3_dv[i];
	}
	dfdx(tmp_x,tmp_v,k4_dx,k4_dv,t+dt,N);

	// Combine to get final dx and dv
	for(size_t i=0U; i<N; ++i){
		dx[i] = one_sixth * (k1_dx[i] + 2.0f * k2_dx[i] + 2.0f * k3_dx[i] + k4_dx[i]);
		dv[i] = one_sixth * (k1_dv[i] + 2.0f * k2_dv[i] + 2.0f * k3_dv[i] + k4_dv[i]);
	}

	// Cleanup
	free(tmp_x); free(tmp_v);
	free(k1_dx); free(k1_dv);
	free(k2_dx); free(k2_dv);
	free(k3_dx); free(k3_dv);
	free(k4_dx); free(k4_dv);
	return;
}

/*
 * Calculates the next 1D coordinates and velocities
 */
void next_1D(float* coord, float* vel, float* new_coord, float* new_vel, float dt, size_t N){
	/* Calculating new coordinates */
	float dx[N];
	float dv[N];
	RK4_1D(coord, vel, dx, dv, 0.0f, dt, &func, N);
	
	for(size_t i=0U; i<N; ++i){
		new_coord[i] = coord[i] + dx[i];
		new_vel[i] = vel[i] + dv[i];
	}
	return;
}


/* --- 2D Structures and Functions ---
                                                                                          
  ▄▄▄▄▄    ▄▄▄▄▄                                                                          
 █▀▀▀▀██▄  ██▀▀▀██                                             ██                         
       ██  ██    ██            ▄▄█████▄  ▀██  ███  ▄▄█████▄  ███████    ▄████▄   ████▄██▄ 
     ▄█▀   ██    ██            ██▄▄▄▄ ▀   ██▄ ██   ██▄▄▄▄ ▀    ██      ██▄▄▄▄██  ██ ██ ██ 
   ▄█▀     ██    ██             ▀▀▀▀██▄    ████▀    ▀▀▀▀██▄    ██      ██▀▀▀▀▀▀  ██ ██ ██ 
 ▄██▄▄▄▄▄  ██▄▄▄██             █▄▄▄▄▄██     ███    █▄▄▄▄▄██    ██▄▄▄   ▀██▄▄▄▄█  ██ ██ ██ 
 ▀▀▀▀▀▀▀▀  ▀▀▀▀▀                ▀▀▀▀▀▀      ██      ▀▀▀▀▀▀      ▀▀▀▀     ▀▀▀▀▀   ▀▀ ▀▀ ▀▀ 
                                          ███                                             
                                                                                          
*/
typedef struct {
	float x;
	float y;
} Vector2D;

void RK4_2D(Vector2D* x, Vector2D* v, Vector2D* dx, Vector2D* dv, float t, float dt,
	    void(*dfdx)(Vector2D*,Vector2D*,Vector2D*,Vector2D*,float,size_t), size_t N){
	/* RK4 Implementation in 2D
	 * x = position array
	 * v = velocity array
	 * dx = derivative of position array
	 * dv = derivative of velocity array
	 * t = current time
	 * dt = time step
	 * dfdx = function that computes derivatives
	 * arguments of dfdx: (x, v, dx, dv, t, N)
	 * N = number of elements
	 */

	// Temporary arrays
	const float one_sixth = 0x1.555556p-3f;
	size_t size = N * sizeof(Vector2D);
	Vector2D* tmp_x = malloc(size);
	Vector2D* tmp_v = malloc(size);

	// k1, k2, k3, k4 arrays for position and velocity
	Vector2D* k1_dx = malloc(size); Vector2D* k1_dv = malloc(size);
	Vector2D* k2_dx = malloc(size); Vector2D* k2_dv = malloc(size);
	Vector2D* k3_dx = malloc(size); Vector2D* k3_dv = malloc(size);
	Vector2D* k4_dx = malloc(size); Vector2D* k4_dv = malloc(size);

	// Calculate k1, k2, k3, k4
	dfdx(x,v,k1_dx,k1_dv,t,N);
	for(size_t i=0U; i<N; ++i){
		tmp_x[i].x = x[i].x + 0.5f * dt * k1_dx[i].x;
		tmp_x[i].y = x[i].y + 0.5f * dt * k1_dx[i].y;
		tmp_v[i].x = v[i].x + 0.5f * dt * k1_dv[i].x;
		tmp_v[i].y = v[i].y + 0.5f * dt * k1_dv[i].y;
	}
	dfdx(tmp_x,tmp_v,k2_dx,k2_dv,t+0.5f*dt,N);
	for(size_t i=0U; i<N; ++i){
		tmp_x[i].x = x[i].x + 0.5f * dt * k2_dx[i].x;
		tmp_x[i].y = x[i].y + 0.5f * dt * k2_dx[i].y;
		tmp_v[i].x = v[i].x + 0.5f * dt * k2_dv[i].x;
		tmp_v[i].y = v[i].y + 0.5f * dt * k2_dv[i].y;
	}
	dfdx(tmp_x,tmp_v,k3_dx,k3_dv,t+0.5f*dt,N);
	for(size_t i=0U; i<N; ++i){
		tmp_x[i].x = x[i].x + dt * k3_dx[i].x;
		tmp_x[i].y = x[i].y + dt * k3_dx[i].y;
		tmp_v[i].x = v[i].x + dt * k3_dv[i].x;
		tmp_v[i].y = v[i].y + dt * k3_dv[i].y;
	}
	dfdx(tmp_x,tmp_v,k4_dx,k4_dv,t+dt,N);

	// Combine to get final dx and dv
	for(size_t i=0U; i<N; ++i){
		dx[i].x = one_sixth * (k1_dx[i].x + 2.0f * k2_dx[i].x + 2.0f * k3_dx[i].x + k4_dx[i].x);
		dx[i].y = one_sixth * (k1_dx[i].y + 2.0f * k2_dx[i].y + 2.0f * k3_dx[i].y + k4_dx[i].y);
		dv[i].x = one_sixth * (k1_dv[i].x + 2.0f * k2_dv[i].x + 2.0f * k3_dv[i].x + k4_dv[i].x);
		dv[i].y = one_sixth * (k1_dv[i].y + 2.0f * k2_dv[i].y + 2.0f * k3_dv[i].y + k4_dv[i].y);
	}

	// Cleanup
	free(tmp_x); free(tmp_v);
	free(k1_dx); free(k1_dv);
	free(k2_dx); free(k2_dv);
	free(k3_dx); free(k3_dv);
	free(k4_dx); free(k4_dv);
	return;
}


/*
 * Calculates the next 2D coordinates and velocities
 */

void next_2D(Vector2D* coord, Vector2D* vel, Vector2D* new_coord, Vector2D* new_vel, float dt, size_t N){
	/* Calculating new coordinates */
	int size = N * sizeof(Vector2D);
	Vector2D* dx = malloc(size);
	Vector2D* dv = malloc(size);
	RK4_2D(coord, vel, dx, dv, 0.0f, dt, &func, N);
	for (size_t i = 0; i < N; i++)
	{
		new_coord[i].x = coord[i].x + dx[i].x * dt;
		new_coord[i].y = coord[i].y + dx[i].y * dt;
		new_vel[i].x = vel[i].x + dv[i].x * dt;
		new_vel[i].y = vel[i].y + dv[i].y * dt;
	}
	free(dx);
	free(dv);
}

/* --- 3D Structures and Functions ---
                                                                                          
  ▄▄▄▄▄    ▄▄▄▄▄                                                                          
 █▀▀▀▀██▄  ██▀▀▀██                                             ██                         
      ▄██  ██    ██            ▄▄█████▄  ▀██  ███  ▄▄█████▄  ███████    ▄████▄   ████▄██▄ 
   █████   ██    ██            ██▄▄▄▄ ▀   ██▄ ██   ██▄▄▄▄ ▀    ██      ██▄▄▄▄██  ██ ██ ██ 
      ▀██  ██    ██             ▀▀▀▀██▄    ████▀    ▀▀▀▀██▄    ██      ██▀▀▀▀▀▀  ██ ██ ██ 
 █▄▄▄▄██▀  ██▄▄▄██             █▄▄▄▄▄██     ███    █▄▄▄▄▄██    ██▄▄▄   ▀██▄▄▄▄█  ██ ██ ██ 
  ▀▀▀▀▀    ▀▀▀▀▀                ▀▀▀▀▀▀      ██      ▀▀▀▀▀▀      ▀▀▀▀     ▀▀▀▀▀   ▀▀ ▀▀ ▀▀ 
                                          ███                                             
*/

typedef struct {
	float x;
	float y;
	float z;
} Vector3D;


void RK4_3D(Vector3D* x, Vector3D* v, Vector3D* dx, Vector3D* dv, float t, float dt,
	    void(*dfdx)(Vector3D*,Vector3D*,Vector3D*,Vector3D*,float,size_t), size_t N){
	/* RK4 Implementation in 3D
	 * x = position array
	 * v = velocity array
	 * dx = derivative of position array
	 * dv = derivative of velocity array
	 * t = current time
	 * dt = time step
	 * dfdx = function that computes derivatives
	 * arguments of dfdx: (x, v, dx, dv, t, N)
	 * N = number of elements
	 */

	// Temporary arrays
	const float one_sixth = 0x1.555556p-3f;
	size_t size = N * sizeof(Vector3D);
	Vector3D* tmp_x = malloc(size);
	Vector3D* tmp_v = malloc(size);

	// k1, k2, k3, k4 arrays for position and velocity
	Vector3D* k1_dx = malloc(size); Vector3D* k1_dv = malloc(size);
	Vector3D* k2_dx = malloc(size); Vector3D* k2_dv = malloc(size);
	Vector3D* k3_dx = malloc(size); Vector3D* k3_dv = malloc(size);
	Vector3D* k4_dx = malloc(size); Vector3D* k4_dv = malloc(size);

	// Calculate k1, k2, k3, k4
	dfdx(x,v,k1_dx,k1_dv,t,N);
	for(size_t i=0U; i<N; ++i){
		tmp_x[i].x = x[i].x + 0.5f * dt * k1_dx[i].x;
		tmp_x[i].y = x[i].y + 0.5f * dt * k1_dx[i].y;
		tmp_x[i].z = x[i].z + 0.5f * dt * k1_dx[i].z;
		tmp_v[i].x = v[i].x + 0.5f * dt * k1_dv[i].x;
		tmp_v[i].y = v[i].y + 0.5f * dt * k1_dv[i].y;
		tmp_v[i].z = v[i].z + 0.5f * dt * k1_dv[i].z;
	}
	dfdx(tmp_x,tmp_v,k2_dx,k2_dv,t+0.5f*dt,N);
	for(size_t i=0U; i<N; ++i){
		tmp_x[i].x = x[i].x + 0.5f * dt * k2_dx[i].x;
		tmp_x[i].y = x[i].y + 0.5f * dt * k2_dx[i].y;
		tmp_x[i].z = x[i].z + 0.5f * dt * k2_dx[i].z;
		tmp_v[i].x = v[i].x + 0.5f * dt * k2_dv[i].x;
		tmp_v[i].y = v[i].y + 0.5f * dt * k2_dv[i].y;
		tmp_v[i].z = v[i].z + 0.5f * dt * k2_dv[i].z;
	}
	dfdx(tmp_x,tmp_v,k3_dx,k3_dv,t+0.5f*dt,N);
	for(size_t i=0U; i<N; ++i){
		tmp_x[i].x = x[i].x + dt * k3_dx[i].x;
		tmp_x[i].y = x[i].y + dt * k3_dx[i].y;
		tmp_x[i].z = x[i].z + dt * k3_dx[i].z;
		tmp_v[i].x = v[i].x + dt * k3_dv[i].x;
		tmp_v[i].y = v[i].y + dt * k3_dv[i].y;
		tmp_v[i].z = v[i].z + dt * k3_dv[i].z;
	}
	dfdx(tmp_x,tmp_v,k4_dx,k4_dv,t+dt,N);

	// Combine to get final dx and dv
	for(size_t i=0U; i<N; ++i){
		dx[i].x = one_sixth * (k1_dx[i].x + 2.0f * k2_dx[i].x + 2.0f * k3_dx[i].x + k4_dx[i].x);
		dx[i].y = one_sixth * (k1_dx[i].y + 2.0f * k2_dx[i].y + 2.0f * k3_dx[i].y + k4_dx[i].y);
		dx[i].z = one_sixth * (k1_dx[i].z + 2.0f * k2_dx[i].z + 2.0f * k3_dx[i].z + k4_dx[i].z);
		dv[i].x = one_sixth * (k1_dv[i].x + 2.0f * k2_dv[i].x + 2.0f * k3_dv[i].x + k4_dv[i].x);
		dv[i].y = one_sixth * (k1_dv[i].y + 2.0f * k2_dv[i].y + 2.0f * k3_dv[i].y + k4_dv[i].y);
		dv[i].z = one_sixth * (k1_dv[i].z + 2.0f * k2_dv[i].z + 2.0f * k3_dv[i].z + k4_dv[i].z);
	}

	// Cleanup
	free(tmp_x); free(tmp_v);
	free(k1_dx); free(k1_dv);
	free(k2_dx); free(k2_dv);
	free(k3_dx); free(k3_dv);
	free(k4_dx); free(k4_dv);
	return;
}


/*
 * Calculates the next 3D coordinates and velocities
 */

void next_3D(Vector3D* coord, Vector3D* vel, Vector3D* new_coord, Vector3D* new_vel, float dt, size_t N){
	/* Calculating new coordinates */
	int size = N * sizeof(Vector3D);
	Vector3D* dx = malloc(size);
	Vector3D* dv = malloc(size);
	RK4_3D(coord, vel, dx, dv, 0.0f, dt, &func, N);
	for (size_t i = 0; i < N; ++i)
	{
		new_coord[i].x = coord[i].x + dx[i].x * dt;
		new_coord[i].y = coord[i].y + dx[i].y * dt;
		new_coord[i].z = coord[i].z + dx[i].z * dt;
		new_vel[i].x = vel[i].x + dv[i].x * dt;
		new_vel[i].y = vel[i].y + dv[i].y * dt;
		new_vel[i].z = vel[i].z + dv[i].z * dt;
	}
	free(dx);
	free(dv);
	return;
}

/*
*
*
	General solver for N generalized coordinates (floats q and dq)
*
*
*/

void RK4_N(float* q, float* dq, float* q_dot, float* dq_dot, float t, float dt,
	    void(*dfdt)(float*,float*,float*,float*,float,size_t), size_t N){
	/* RK4 Implementation in N dimensions
	 * q = generalized coordinates array
	 * dq = generalized velocities array
	 * dq_out = derivative of generalized coordinates array
	 * ddq_out = derivative of generalized velocities array
	 * t = current time
	 * dt = time step
	 * dfdt = function that computes derivatives
	 * arguments of dfdt: (q, dq, _dq, _ddq, t, N)
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
	// Combine to get final dq and ddq
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
}