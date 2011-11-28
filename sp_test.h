#ifndef _SP_TEST_H_
#define _SP_TEST_H_

#define TEST(fname) static int test_#fname()
#define ADD_TEST(name) printf("test_name:\t*%s*\n",(test_#test_name) ? "pass" : "fail"); 

#endif /* _SP_TEST_H_ */
