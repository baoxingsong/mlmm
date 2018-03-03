//
// Created by baoxing on 2/26/18.
//

#include "../../impl/my_math_impl.h"
#include "../../../googletest/googletest/include/gtest/gtest.h"


TEST(sf, c1){
    int n1,n2,i;
    double y,f;
    static int n[2]={ 2,5};
    static int m[2]={ 3,10};
    printf("\n");
    for (i=0; i<=1; i++)
    { n1=n[i]; n2=m[i]; f=3.5;
        y=sf(f,n1,n2);
        printf("P(%4.2f, %d, %d)=%e\n",f,n1,n2,y);
        f=9.0; y=sf(f,n1,n2);
        printf("P(%4.2f, %d, %d)=%e\n",f,n1,n2,y);
    }
    printf("\n");

    f = 10;
    n1 = 2;
    n2 = 40;
    std::cout << sf(f, n1, n2) << std::endl;
    ASSERT_EQ(0, 0);
}
