/*
   

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
*/

#pragma once

#include <stddef.h>

int calloc_s(void **ptr, size_t nelem, size_t elsize);

int malloc_s(void **ptr, size_t size);

int realloc_s(void **ptr, size_t new_size);