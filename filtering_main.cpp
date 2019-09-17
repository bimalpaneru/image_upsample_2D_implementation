/*****************************************************************************/
// File: filtering_main.cpp
// Author: David Taubman
// Last Revised: 13 August, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

#include "io_bmp.h"
#include "image_comps.h"
#include<math.h>
#include<cmath>
#include<time.h>

const float PI = 3.14F;




/* ========================================================================= */
/*                 Implementation of `my_image_comp' functions               */
/* ========================================================================= */

/*****************************************************************************/
/*                  my_image_comp::perform_boundary_extension                */
/*****************************************************************************/

float sinc(float sample) {
	if (sample == 0)
		return 1.0F;
	else
		return (sinf(PI * sample ) / (PI * sample));

}

float hann_win(float x, int wnd_size) {
	//float sum = 0.5F * (1.0F + cosf(PI * x / float(wnd_size)));
	return 0.5F * (1.0F + cosf(PI * x / float(wnd_size)));
}


void matrix_mul(float* mat_1, float* mat_2, float* result, int filter_H) {
	float sum = 0.0F;
	for (int x = 0; x <= 2 * filter_H; x++) {
		for (int y = 0; y <= 2 * filter_H; y++) {
			result[x * (2 * filter_H + 1) + y] = mat_1[x] * mat_2[y];
			sum += (result[x * (2 * filter_H + 1) + y]);
			
		}

	}
	//printf("gain= %f\n", sum);
	for (int i = 0; i < (2 * filter_H + 1)*(2 * filter_H + 1); i++) {

		result[i] = result[i] / sum;
	}

}

//if transpose == 1 then transpose filter matrix if 0 then don't
float inner_product(float* ip, int ip_stride, float* mirror_psf, int filter_extent, int transpose) {
	float sum = 0.0F;
	int x_t = 0;
	int y_t = 1;
	if (transpose == 1) {
		y_t = 0;
		x_t = 1;

	}
	else {
		y_t = 1;
		x_t = 0;

	}
	//printf("Transpose matrix \n");
	for (int y = -filter_extent; y <= filter_extent; y++) {
		for (int x = -filter_extent; x <= filter_extent; x++) {
			//ip[y * ip_stride + x];
			//mirror_psf[y * (2 * filter_extent + 1) + x];
			//printf(" %f \t ", mirror_psf[y * (2 * filter_extent * y_t + 1) + x * (2 * filter_extent * x_t + 1)]);
			sum += ip[y * ip_stride + x] * mirror_psf[y * (2 * filter_extent * y_t + 1) + x * (2 * filter_extent * x_t + 1) ];
			
		}
		//printf("\n");
	}
	return sum;
}



void my_image_comp::perform_boundary_extension()
{
  int r, c;

  // First extend upwards
  float *first_line = buf;
  for (r = 1; r <= border; r++)
	  for (c = 0; c < width; c++)
		  first_line[-r * stride + c] = first_line[(r - 1) * stride + c];//first_line[c];

  // Now extend downwards
  float *last_line = buf+(height-1)*stride;
  for (r = 1; r <= border; r++)
	  for (c = 0; c < width; c++)
		  last_line[r * stride + c] = last_line[-(r - 1) * stride + c];//[(height - 1) * -r * stride + c];    //     r * stride + c];//0;//last_line[c];

  // Now extend all rows to the left and to the right
  float *left_edge = buf-border*stride;
  float *right_edge = left_edge + width - 1;
  for (r=height+2*border; r > 0; r--, left_edge+=stride, right_edge+=stride)
    for (c=1; c <= border; c++)
      {
		left_edge[-c] =  left_edge[c];// left_edge[0];
		right_edge[c] =  right_edge[-c];// right_edge[0];
      }
}


/* ========================================================================= */
/*                              Global Functions                             */
/* ========================================================================= */

/*****************************************************************************/
/*                                apply_filter                               */
/*****************************************************************************/

void get_filter_kernels(float* q_5m_p0,float* q_5m_p1, float* q_5m_p2, float* q_5m_p3, float* q_5m_p4, int filter_H) {
	/*float s0 = 0.0F;
	float s1 = 0.0F;
	float s2 = 0.0F;
	float s3 = 0.0F;
	float s4 = 0.0F;*/

	for (int i = 0 ; i <= 2*filter_H; i++) {
		
		q_5m_p0[i] = sinc(float(i) - float(filter_H)) *   hann_win(float(i) - float(filter_H), filter_H + 1);
		q_5m_p1[i] = sinc(float(i) - float(filter_H) - 0.4F) * hann_win(float(i) - float(filter_H) - 0.4F, filter_H + 1 );
		q_5m_p2[i] = sinc(float(i) - float(filter_H) + 0.2F) * hann_win(float(i) - float(filter_H) + 0.2F, filter_H + 1);
		q_5m_p3[i] = sinc(float(i) - float(filter_H) - 0.2F) * hann_win(float(i) - float(filter_H) - 0.2F, filter_H + 1);
		q_5m_p4[i] = sinc(float(i) - float(filter_H) + 0.4F) * hann_win(float(i) - float(filter_H) + 0.4F, filter_H + 1);
		/*s0 += q_5m_p0[i];
		s1 += q_5m_p1[i];
		s2 += q_5m_p2[i];
		s3 += q_5m_p3[i];
		s4 += q_5m_p4[i];*/
	    //printf("%f \t", q_5m_p0[i]);
	
	}

	/*for (int i = 0; i <= 2 * filter_H; i++) {

		q_5m_p0[i] = q_5m_p0[i] / s0;
		q_5m_p1[i] = q_5m_p1[i] / s1;
		q_5m_p2[i] = q_5m_p2[i] / s2;
		q_5m_p3[i] = q_5m_p3[i] / s3;
		q_5m_p4[i] = q_5m_p4[i] / s4;
		
		//printf("%f \t", q_5m_p0[i]);

	}*/
	//for (int i = 0; i <= 2 * filter_H; i++)
	//printf("%f \t", q_5m_p0[i]);




       // printf("\n");

  //for (int i = 0; i <= 2 * filter_H; i++)
    //   printf("%f \t", q_5m_p1[i]);

}
void apply_filter(my_image_comp *in, my_image_comp *out, int filter_H)
{
	   
	 //filter length H in 1D then 2H+1 is 1D filter length
	float* q_5m_p0 = new float[2 * filter_H + 1];
	float* q_5m_p1 = new float[2 * filter_H + 1];
	float* q_5m_p2 = new float[2 * filter_H + 1];
	float* q_5m_p3 = new float[2 * filter_H + 1];
	float* q_5m_p4 = new float[2 * filter_H + 1];

	get_filter_kernels(q_5m_p0,q_5m_p1, q_5m_p2, q_5m_p3, q_5m_p4, filter_H);

	float* m5_p0_5m_p1 = new float[(2 * filter_H + 1) * (2 * filter_H + 1)];
	float* m5_p1_5m_p1 = new float[(2 * filter_H + 1) * (2 * filter_H + 1)];
	float* m5_p1_5m_p2 = new float[(2 * filter_H + 1) * (2 * filter_H + 1)];
	float* m5_p1_5m_p3 = new float[(2 * filter_H + 1) * (2 * filter_H + 1)];
	float* m5_p1_5m_p4 = new float[(2 * filter_H + 1) * (2 * filter_H + 1)];
	
	float* m5_p0_5m_p2 = new float[(2 * filter_H + 1) * (2 * filter_H + 1)];
	float* m5_p2_5m_p2 = new float[(2 * filter_H + 1) * (2 * filter_H + 1)];
	float* m5_p2_5m_p3 = new float[(2 * filter_H + 1) * (2 * filter_H + 1)];
	float* m5_p2_5m_p4 = new float[(2 * filter_H + 1) * (2 * filter_H + 1)];
	
	float* m5_p0_5m_p3 = new float[(2 * filter_H + 1) * (2 * filter_H + 1)];
	float* m5_p3_5m_p3 = new float[(2 * filter_H + 1) * (2 * filter_H + 1)];
	float* m5_p3_5m_p4 = new float[(2 * filter_H + 1) * (2 * filter_H + 1)];
	
	float* m5_p0_5m_p4 = new float[(2 * filter_H + 1) * (2 * filter_H + 1)];
	float* m5_p4_5m_p4 = new float[(2 * filter_H + 1) * (2 * filter_H + 1)];
	
	matrix_mul(q_5m_p0, q_5m_p1, m5_p0_5m_p1, filter_H);
	matrix_mul(q_5m_p1, q_5m_p1, m5_p1_5m_p1, filter_H);
	matrix_mul(q_5m_p1, q_5m_p2, m5_p1_5m_p2, filter_H);
	matrix_mul(q_5m_p1, q_5m_p3, m5_p1_5m_p3, filter_H);
	matrix_mul(q_5m_p1, q_5m_p4, m5_p1_5m_p4, filter_H);

	matrix_mul(q_5m_p0, q_5m_p2, m5_p0_5m_p2, filter_H);
	matrix_mul(q_5m_p2, q_5m_p2, m5_p2_5m_p2, filter_H);
	matrix_mul(q_5m_p2, q_5m_p3, m5_p2_5m_p3, filter_H);
	matrix_mul(q_5m_p2, q_5m_p4, m5_p2_5m_p4, filter_H);

	matrix_mul(q_5m_p0, q_5m_p3, m5_p0_5m_p3, filter_H);
	matrix_mul(q_5m_p3, q_5m_p3, m5_p3_5m_p3, filter_H);
	matrix_mul(q_5m_p3, q_5m_p4, m5_p3_5m_p4, filter_H);

	matrix_mul(q_5m_p0, q_5m_p4, m5_p0_5m_p4, filter_H);
	matrix_mul(q_5m_p4, q_5m_p4, m5_p4_5m_p4, filter_H);

	/*for (int x = 0; x <= 2 * filter_H; x++) {

		for (int y = 0; y <= 2 * filter_H; y++) {
			printf("%f \t", m5_p0_5m_p1[x * (2 * filter_H + 1) + y]);

		}
		printf("\n");

	}*/

	


	
	m5_p0_5m_p1 = m5_p0_5m_p1 + (2 * filter_H + 1) * filter_H + filter_H;
	m5_p1_5m_p1 = m5_p1_5m_p1 + (2 * filter_H + 1) * filter_H + filter_H;
	m5_p1_5m_p2 = m5_p1_5m_p2 + (2 * filter_H + 1) * filter_H + filter_H;
    m5_p1_5m_p3 = m5_p1_5m_p3 + (2 * filter_H + 1) * filter_H + filter_H;
    m5_p1_5m_p4 = m5_p1_5m_p4 + (2 * filter_H + 1) * filter_H + filter_H;

	m5_p0_5m_p2 = m5_p0_5m_p2 + (2 * filter_H + 1) * filter_H + filter_H;
	m5_p2_5m_p2 = m5_p2_5m_p2 + (2 * filter_H + 1) * filter_H + filter_H;
	m5_p2_5m_p3 = m5_p2_5m_p3 + (2 * filter_H + 1) * filter_H + filter_H;
	m5_p2_5m_p4 = m5_p2_5m_p4 + (2 * filter_H + 1) * filter_H + filter_H;

	m5_p0_5m_p3 = m5_p0_5m_p3 + (2 * filter_H + 1) * filter_H + filter_H;
	m5_p3_5m_p3 = m5_p3_5m_p3 + (2 * filter_H + 1) * filter_H + filter_H;
	m5_p3_5m_p4 = m5_p3_5m_p4 + (2 * filter_H + 1) * filter_H + filter_H;

	m5_p0_5m_p4 = m5_p0_5m_p4 + (2 * filter_H + 1) * filter_H + filter_H;
	m5_p4_5m_p4 = m5_p4_5m_p4 + (2 * filter_H + 1) * filter_H + filter_H;

	
	//float* ip = in->buf + 1 * in->stride + 3;
	//inner_product(ip, in->stride, m5_p1_5m_p2, filter_H, 1);




	    
  // Check for consistent dimensions
  assert(in->border >= filter_H);
  //assert((out->height <= new_width) && (out->width <= in->width));
  
  // Perform the convolution
  for (int r = 0; r < out->height; r++) {

	  for (int c = 0; c < out->width; c++) {

		  float* op = out->buf + r * out->stride + c;


		  //r mod 5 == 0
		  if (r % 5 == 0 && c % 5 == 0) {
			  int r_ip = r * 2 / 5;
			  int c_ip = c * 2 / 5;
			  float* ip = in->buf + r_ip * in->stride + c_ip;
			  *op = *ip;

		  }
		 
		   else if (r % 5 == 0 && c % 5 == 1) {
			  int r_ip = r * 2 / 5;
			  int c_ip = (c - 1) * 2 / 5;
			  float* ip = in->buf + r_ip * in->stride + c_ip;
			  *op = inner_product(ip, in->stride, m5_p0_5m_p1, filter_H, 0);

		  }
		   
		   

		  
		  else if (r % 5 == 0 && c % 5 == 2) {
			  int r_ip = r * 2 / 5;
			  int c_ip = (c - 2) * 2 / 5 + 1;
			  float* ip = in->buf + r_ip * in->stride + c_ip;
			  *op = inner_product(ip, in->stride, m5_p0_5m_p2, filter_H, 0);

		  }


		  
			
		  else if (r % 5 == 0 && c % 5 == 3) {
			  int r_ip = r * 2 / 5;
			  int c_ip = (c - 3) * 2 / 5 + 1;
			  float* ip = in->buf + r_ip * in->stride + c_ip;
			  *op = inner_product(ip, in->stride, m5_p0_5m_p3, filter_H, 0);

		  }
		  
		 
		 else if (r % 5 == 0 && c % 5 == 4) {
			  int r_ip = r * 2 / 5;
			  int c_ip = (c - 4) * 2 / 5 + 2;
			  float* ip = in->buf + r_ip * in->stride + c_ip;
			  *op = inner_product(ip, in->stride, m5_p0_5m_p4, filter_H, 0);

		  }


		  // r mod 5 = 1 from here
		  else if (r % 5 == 1 && c % 5 == 0) {
			  int r_ip = (r - 1) * 2 / 5;
			  int c_ip = c * 2 / 5;
			  float* ip = in->buf + r_ip * in->stride + c_ip;
			  *op = inner_product(ip, in->stride, m5_p0_5m_p1, filter_H, 1);

		  }

		  else if (r % 5 == 1 && c % 5 == 1) {
			  int r_ip = (r - 1) * 2 / 5;
			  int c_ip = (c - 1) * 2 / 5;
			  float* ip = in->buf + r_ip * in->stride + c_ip;
			  *op = inner_product(ip, in->stride, m5_p1_5m_p1, filter_H, 0);

		  }

		  else if (r % 5 == 1 && c % 5 == 2) {
			  int r_ip = (r - 1) * 2 / 5;
			  int c_ip = (c - 2) * 2 / 5 + 1;
			  float* ip = in->buf + r_ip * in->stride + c_ip;
			  *op = inner_product(ip, in->stride, m5_p1_5m_p2, filter_H, 0);

		  }

		   else if (r % 5 == 1 && c % 5 == 3) {
			  int r_ip = (r - 1) * 2 / 5;
			  int c_ip = (c - 3) * 2 / 5 + 1;
			  float* ip = in->buf + r_ip * in->stride + c_ip;
			  *op = inner_product(ip, in->stride, m5_p1_5m_p3, filter_H, 0);

		  }

		 	  
		  else if (r % 5 == 1 && c % 5 == 4) {
			  int r_ip = (r - 1) * 2 / 5;
			  int c_ip = (c - 4) * 2 / 5 + 2;
			  float* ip = in->buf + r_ip * in->stride + c_ip;
			  *op = inner_product(ip, in->stride, m5_p1_5m_p4, filter_H, 0);

		  }


		  //rmod 5 == 2
		  else if (r % 5 == 2 && c % 5 == 0) {
			  int r_ip = (r - 2) * 2 / 5 + 1;
			  int c_ip = c * 2 / 5;
			  float* ip = in->buf + r_ip * in->stride + c_ip;
			  *op = inner_product(ip, in->stride, m5_p0_5m_p2, filter_H, 1);

		  }

		  else if (r % 5 == 2 && c % 5 == 1) {
			  int r_ip = (r - 2) * 2 / 5 + 1;
			  int c_ip = (c - 1) * 2 / 5;
			  float* ip = in->buf + r_ip * in->stride + c_ip;
			  *op = inner_product(ip, in->stride, m5_p1_5m_p2, filter_H, 1); //p2-p1 by transpose

		  }

		  else if (r % 5 == 2 && c % 5 == 2) {
			  int r_ip = (r - 2) * 2 / 5 + 1;
			  int c_ip = (c - 2) * 2 / 5 + 1;
			  float* ip = in->buf + r_ip * in->stride + c_ip;
			  *op = inner_product(ip, in->stride, m5_p2_5m_p2, filter_H, 0);

		  }

		  else if (r % 5 == 2 && c % 5 == 3) {
			  int r_ip = (r - 2) * 2 / 5 + 1;
			  int c_ip = (c - 3) * 2 / 5 + 1;
			  float* ip = in->buf + r_ip * in->stride + c_ip;
			  *op = inner_product(ip, in->stride, m5_p2_5m_p3, filter_H, 0);

		  }


		  else if (r % 5 == 2 && c % 5 == 4) {
			  int r_ip = (r - 2) * 2 / 5 + 1;
			  int c_ip = (c - 4) * 2 / 5 + 2;
			  float* ip = in->buf + r_ip * in->stride + c_ip;
			  *op = inner_product(ip, in->stride, m5_p2_5m_p4, filter_H, 0);

		  }

		  //rmod 5 == 3
		  else if (r % 5 == 3 && c % 5 == 0) {
			  int r_ip = (r - 3) * 2 / 5 + 1;
			  int c_ip = c * 2 / 5;
			  float* ip = in->buf + r_ip * in->stride + c_ip;
			  *op = inner_product(ip, in->stride, m5_p0_5m_p3, filter_H, 1);

		  }

		  else if (r % 5 == 3 && c % 5 == 1) {
			  int r_ip = (r - 3) * 2 / 5 + 1;
			  int c_ip = (c - 1) * 2 / 5;
			  float* ip = in->buf + r_ip * in->stride + c_ip;
			  *op = inner_product(ip, in->stride, m5_p1_5m_p3, filter_H, 1); //p3-p1 by transpose

		  }

		  else if (r % 5 == 3 && c % 5 == 2) {
			  int r_ip = (r - 3) * 2 / 5 + 1;
			  int c_ip = (c - 2) * 2 / 5 + 1;
			  float* ip = in->buf + r_ip * in->stride + c_ip;
			  *op = inner_product(ip, in->stride, m5_p2_5m_p3, filter_H, 1);

		  }

		  else if (r % 5 == 3 && c % 5 == 3) {
			  int r_ip = (r - 3) * 2 / 5 + 1;
			  int c_ip = (c - 3) * 2 / 5 + 1;
			  float* ip = in->buf + r_ip * in->stride + c_ip;
			  *op = inner_product(ip, in->stride, m5_p3_5m_p3, filter_H, 0);

		  }


		  else if (r % 5 == 3 && c % 5 == 4) {
			  int r_ip = (r - 3) * 2 / 5 + 1;
			  int c_ip = (c - 4) * 2 / 5 + 2;
			  float* ip = in->buf + r_ip * in->stride + c_ip;
			  *op = inner_product(ip, in->stride, m5_p3_5m_p4, filter_H, 0);

		  }


		  //rmod 5 == 4
		  else if (r % 5 == 4 && c % 5 == 0) {
			  int r_ip = (r - 4) * 2 / 5 + 2;
			  int c_ip = c * 2 / 5;
			  float* ip = in->buf + r_ip * in->stride + c_ip;
			  *op = inner_product(ip, in->stride, m5_p0_5m_p4, filter_H, 1);

		  }

		  else if (r % 5 == 4 && c % 5 == 1) {
			  int r_ip = (r - 4) * 2 / 5 + 2;
			  int c_ip = (c - 1) * 2 / 5;
			  float* ip = in->buf + r_ip * in->stride + c_ip;
			  *op = inner_product(ip, in->stride, m5_p1_5m_p4, filter_H, 1); //p3-p1 by transpose

		  }

		  else if (r % 5 == 4 && c % 5 == 2) {
			  int r_ip = (r - 4) * 2 / 5 + 2;
			  int c_ip = (c - 2) * 2 / 5 + 1;
			  float* ip = in->buf + r_ip * in->stride + c_ip;
			  *op = inner_product(ip, in->stride, m5_p2_5m_p4, filter_H, 1);

		  }

		  else if (r % 5 == 4 && c % 5 == 3) {
			  int r_ip = (r - 4) * 2 / 5 + 2;
			  int c_ip = (c - 3) * 2 / 5 + 1;
			  float* ip = in->buf + r_ip * in->stride + c_ip;
			  *op = inner_product(ip, in->stride, m5_p3_5m_p4, filter_H, 1);

		  }


		  else if (r % 5 == 4 && c % 5 == 4) {
			  int r_ip = (r - 4) * 2 / 5 + 2;
			  int c_ip = (c - 4) * 2 / 5 + 2;
			  float* ip = in->buf + r_ip * in->stride + c_ip;
			  *op = inner_product(ip, in->stride, m5_p4_5m_p4, filter_H, 0);

		  }
		

	  }
  }


   
}

/*****************************************************************************/
/*                                    main                                   */
/*****************************************************************************/

int  main(int argc, char *argv[])
{
	clock_t start_time = clock();
 
   
  if (argc != 4)
    {
      fprintf(stderr,"Usage: %s <in bmp file> <out bmp file>\n",argv[0]);
      return -1;
    }
  int filter_H = atoi(argv[3]);
 // printf("filter length %d \n", filter_H);
 
  int err_code=0;
  try {
      // Read the input image
      bmp_in in;
      if ((err_code = bmp_in__open(&in,argv[1])) != 0)
        throw err_code;

	  
	
      int width = in.cols, height = in.rows;

	  int new_width = int(ceil(float(width) * 2.5F));
	  int new_height = int(ceil(float(height) *2.5F));

      int n, num_comps = in.num_components;
      my_image_comp *input_comps = new my_image_comp[num_comps];
      
	  ///extent boundaries
	  for (n = 0; n < num_comps; n++) {
		  input_comps[n].init(height, width, filter_H); // Leave a border of filter extent size
	  }
        
      

      int r; // Declare row index
      io_byte *line = new io_byte[width*num_comps];
      for (r=height-1; r >= 0; r--)
        { // "r" holds the true row index we are reading, since the image is
          // stored upside down in the BMP file.
          if ((err_code = bmp_in__get_line(&in,line)) != 0)
            throw err_code;
          for (n=0; n < num_comps; n++)
            {
              io_byte *src = line+n; // Points to first sample of component n
              float *dst = input_comps[n].buf + r * input_comps[n].stride;
              for (int c=0; c < width; c++, src+=num_comps)
            				  
				  dst[c] = (float) *src; // The cast to type "float" is not
                      // strictly required here, since bytes can always be
                      // converted to floats without any loss of information.
            }
        }
      bmp_in__close(&in);

      // Allocate storage for the filtered output
      my_image_comp *output_comps = new my_image_comp[num_comps];

	  for (n = 0; n < num_comps; n++) {
		  output_comps[n].init(new_height, new_width, 0); // Don't need a border for output
	  }
	  

      // Process the image, all in floating point (easy)
     for (n=0; n < num_comps; n++)
        input_comps[n].perform_boundary_extension();
      for (n=0; n < num_comps; n++)
        apply_filter(input_comps+n,output_comps+n, filter_H);

      // Write the image back out again
     bmp_out out;
	 io_byte* lineout = new io_byte[new_width * num_comps];
     if ((err_code = bmp_out__open(&out,argv[2],new_width,new_height,num_comps)) != 0)
        throw err_code;
      for (r=new_height-1; r >= 0; r--)
        { // "r" holds the true row index we are writing, since the image is
          // written upside down in BMP files.
          for (n=0; n < num_comps; n++)
            {
              io_byte *dst = lineout+n; // Points to first sample of component n
              float *src = output_comps[n].buf + r * output_comps[n].stride;
			  for (int c = 0; c < new_width; c++, dst += num_comps)
			  {
				  if (src[c] > 255.0F)
					  src[c] = 255;
				  else if (src[c] <0.0F)
					  src[c] = 0;

				  *dst = (io_byte)src[c];
			  } // The cast to type "io_byte" is
                      // required here, since floats cannot generally be
                      // converted to bytes without loss of information.  The
                      // compiler will warn you of this if you remove the cast.
                      // There is in fact not the best way to do the
                      // conversion.  You should fix it up in the lab.
            }
          bmp_out__put_line(&out,lineout);
        }
      bmp_out__close(&out);
      
      delete[] input_comps;
      delete[] output_comps;
    }
  catch (int exc) {
      if (exc == IO_ERR_NO_FILE)
        fprintf(stderr,"Cannot open supplied input or output file.\n");
      else if (exc == IO_ERR_FILE_HEADER)
        fprintf(stderr,"Error encountered while parsing BMP file header.\n");
      else if (exc == IO_ERR_UNSUPPORTED)
        fprintf(stderr,"Input uses an unsupported BMP file format.\n  Current "
                "simple example supports only 8-bit and 24-bit data.\n");
      else if (exc == IO_ERR_FILE_TRUNC)
        fprintf(stderr,"Input or output file truncated unexpectedly.\n");
      else if (exc == IO_ERR_FILE_NOT_OPEN)
        fprintf(stderr,"Trying to access a file which is not open!(?)\n");
      return -1;
    }
  clock_t end_time = clock();
  float seconds = (float((end_time - start_time))) / CLOCKS_PER_SEC;
  printf("Computation Time = %f\n", seconds);
  return 0;
}
