//
//  checks.h
//  
//
//  Created by Johannes Dietschreit on 11.06.15.
//
//

#ifndef ____checks__
#define ____checks__

int check_saw1(int (**check_ptr), int position);

int check_saw2(int (**check_ptr), int position, int lastmove_var, const int (sawmoves_ptr)[4][3] , const int (moves_ptr)[2][4][3], int odd_var );

int check_saw3(int (**check_ptr), int position, int lastmove_var, const int (sawmoves_ptr)[4][3] , const int (moves_ptr)[2][4][3], int odd_var );

int check_fsaw3(int (**check_ptr), int position, int start_f, int lastmove_var, const int (sawmoves_ptr)[4][3] , const int (moves_ptr)[2][4][3], int odd_var );

int check_whole_fsaw3(int (**check_ptr), int start_f, int end_f, const int (moves_ptr)[2][4][3]);

int check_bsaw3(int (**check_ptr), int position, int lastmove_var, const int (sawmoves_ptr)[4][3] , const int (moves_ptr)[2][4][3], int odd_var , int bpos_var);

int check_fb_saw3(int (**check_ptr), int position, int start_f, int lastmove_var, const int (sawmoves_ptr)[4][3] , const int (moves_ptr)[2][4][3], int odd_var , int bpos_var);


#endif /* defined(____checks__) */
