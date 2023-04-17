#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

#include "test.h"

int main()
{
    CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( PolyTest , "t1" );
    CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( InterpTest , "t2" );
    CppUnit::TextUi::TestRunner runner1;
    CppUnit::TextUi::TestRunner runner2;
    CppUnit::TestFactoryRegistry &registry1 = CppUnit::TestFactoryRegistry::getRegistry("t1");
    CppUnit::TestFactoryRegistry &registry2 = CppUnit::TestFactoryRegistry::getRegistry("t2");
    runner1.addTest(registry1.makeTest());
    runner2.addTest(registry2.makeTest());
    runner1.run();
    runner2.run();
}
