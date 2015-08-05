/*
   A monotonic clock which increases by 1 second for every real second, can be
   used to time program execution.

   Copyright 2015 Ian Kavanagh

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   Adapted from:
   ThomasHabets/monotonic_clock <thomas@habets.se> 2010
   https://github.com/ThomasHabets/monotonic_clock
*/

#pragma once

/**
 * Purpose
 * =======
 *
 * To time the execution of segments of code. Uses a monotonic clock which
 * increases by 1 for every 1 second passed in time.
 * 
 * Output
 * ======
 *
 *          DOUBLE.
 *          A number of seconds.
 *
 * =============================================================================
 */
double monotonic_clock(void);
