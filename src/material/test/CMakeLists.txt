#------------------------------------------------------------------------------#
# libdetran/src/material/test/CMakeLists.txt
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Create test executables
#------------------------------------------------------------------------------#

ADD_EXECUTABLE(test_Material test_Material.cc)
ADD_EXECUTABLE(test_dummy_material test_dummy_material.cc)

TARGET_LINK_LIBRARIES(test_Material material utilities)
TARGET_LINK_LIBRARIES(test_dummy_material material utilities)

#------------------------------------------------------------------------------#
# Define tests
#------------------------------------------------------------------------------#

ADD_TEST( test_Material_basic  test_Material 0)
ADD_TEST( test_Material_bounds test_Material 1)
#ADD_TEST(test_dummy_material test_Material 3)
ADD_TEST(test_dummy_material test_dummy_material  0)


