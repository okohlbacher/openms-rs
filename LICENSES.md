# Licenses for OpenMS for Rust

The OpenMS-derived implementation is BSD-3-Clause, as reproduced below.
The private random helper shared by decoy generation and unique IDs additionally retains the Boost Software License
1.0 and its original author notices; its component terms are reproduced below.
The embedded OpenMS Rust Modification Table includes transformed UniMod
data under the Design Science License. The original source XML, transformation
script, data notices and full license are supplied in resources/modifications.
These are component-specific terms, not a choice of license for the same work.

The bundled XLMOD 2016 vocabulary and its transformed reference table are
CC BY 3.0, attributed to Lutz Fischer and Gerhard Mayer. Its version, source,
license link and unchanged-data notice are preserved in
[the XLMOD notice](resources/modifications/XLMOD_LICENSE.md).

The pinned OpenMS enzyme registry in `resources/enzymes/Enzymes.xml` and its
generated Rust predicates are included under the OpenMS BSD-3-Clause terms.

The RNA snapshot and its embedded projection retain separate MODOMICS data
notices in [the RNA data record](resources/rna/README.md). General redistribution
terms for this dataset are unresolved; the software BSD license does not grant
rights to the MODOMICS data. The complete snapshot is retained for the port
and its scientific validation; public source availability does not resolve
the separate dataset redistribution terms.

External Cargo dependencies retain their own licenses and notices. The EMG
fitter uses `libm` 0.2.16 (MIT) for the complementary error function; its source
and full notices are distributed in the dependency's Cargo package. No `libm`
source is vendored in this repository. The optional MODOMICS JSON reader uses
serde_json 1.0.150 (MIT OR Apache-2.0), with its own Cargo dependency notices.
Optional file compression uses bzip2 0.6.1 (MIT OR Apache-2.0) with its default
pure Rust libbz2-rs-sys backend (bzip2-1.0.6), and flate2 1.1.9
(MIT OR Apache-2.0). Their full licenses are supplied in the Cargo packages;
no compression-library source is vendored here.

## OpenMS implementation and custom data: BSD-3-Clause

--------------------------------------------------------------------------
                  OpenMS -- Open-Source Mass Spectrometry
--------------------------------------------------------------------------
Copyright OpenMS Inc. -- Eberhard Karls University Tuebingen,
ETH Zurich, and Freie Universitaet Berlin 2002-present.

This software is released under a three-clause BSD license:
 * Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
 * Neither the name of any author or any participating institution
   may be used to endorse or promote products derived from this software
   without specific prior written permission.
For a full list of authors, refer to the file AUTHORS.
--------------------------------------------------------------------------
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


## UniMod data and its derivative: Design Science License

Copyright (C) 2002-2006 Unimod.

DESIGN SCIENCE LICENSE

TERMS AND CONDITIONS FOR COPYING, DISTRIBUTION AND MODIFICATION

Copyright © 1999-2001 Michael Stutz <stutz@dsl.org>
Verbatim copying of this document is permitted, in any medium.

0. PREAMBLE.

Copyright law gives certain exclusive rights to the author of a work,
including the rights to copy, modify and distribute the work (the
"reproductive," "adaptative," and "distribution" rights).

The idea of "copyleft" is to willfully revoke the exclusivity of those
rights under certain terms and conditions, so that anyone can copy and
distribute the work or properly attributed derivative works, while all
copies remain under the same terms and conditions as the original.

The intent of this license is to be a general "copyleft" that can be
applied to any kind of work that has protection under copyright. This
license states those certain conditions under which a work published
under its terms may be copied, distributed, and modified.

Whereas "design science" is a strategy for the development of
artifacts as a way to reform the environment (not people) and
subsequently improve the universal standard of living, this Design
Science License was written and deployed as a strategy for promoting
the progress of science and art through reform of the environment.

1. DEFINITIONS.

"License" shall mean this Design Science License. The License applies
to any work which contains a notice placed by the work's copyright
holder stating that it is published under the terms of this Design
Science License.

"Work" shall mean such an aforementioned work. The License also
applies to the output of the Work, only if said output constitutes a
"derivative work" of the licensed Work as defined by copyright law.

"Object Form" shall mean an executable or performable form of the
Work, being an embodiment of the Work in some tangible medium.

"Source Data" shall mean the origin of the Object Form, being the
entire, machine-readable, preferred form of the Work for copying and
for human modification (usually the language, encoding or format in
which composed or recorded by the Author); plus any accompanying
files, scripts or other data necessary for installation, configuration
or compilation of the Work.

(Examples of "Source Data" include, but are not limited to, the
following: if the Work is an image file composed and edited in PNG
format, then the original PNG source file is the Source Data; if the
Work is an MPEG 1.0 layer 3 digital audio recording made from a WAV
format audio file recording of an analog source, then the original WAV
file is the Source Data; if the Work was composed as an unformatted
plaintext file, then that file is the Source Data; if the Work was
composed in LaTeX, the LaTeX file(s) and any image files and/or custom
macros necessary for compilation constitute the Source Data.)

"Author" shall mean the copyright holder(s) of the Work.

The individual licensees are referred to as "you."

2. RIGHTS AND COPYRIGHT.

The Work is copyrighted by the Author. All rights to the Work are
reserved by the Author, except as specifically described below. This
License describes the terms and conditions under which the Author
permits you to copy, distribute and modify copies of the Work.

In addition, you may refer to the Work, talk about it, and (as
dictated by "fair use") quote from it, just as you would any
copyrighted material under copyright law.

Your right to operate, perform, read or otherwise interpret and/or
execute the Work is unrestricted; however, you do so at your own risk,
because the Work comes WITHOUT ANY WARRANTY -- see Section 7 ("NO
WARRANTY") below.

3. COPYING AND DISTRIBUTION.

Permission is granted to distribute, publish or otherwise present
verbatim copies of the entire Source Data of the Work, in any medium,
provided that full copyright notice and disclaimer of warranty, where
applicable, is conspicuously published on all copies, and a copy of
this License is distributed along with the Work.

Permission is granted to distribute, publish or otherwise present
copies of the Object Form of the Work, in any medium, under the terms
for distribution of Source Data above and also provided that one of
the following additional conditions are met:

(a) The Source Data is included in the same distribution, distributed
under the terms of this License; or

(b) A written offer is included with the distribution, valid for at
least three years or for as long as the distribution is in print
(whichever is longer), with a publicly-accessible address (such as a
URL on the Internet) where, for a charge not greater than
transportation and media costs, anyone may receive a copy of the
Source Data of the Work distributed according to the section above; or

(c) A third party's written offer for obtaining the Source Data at no
cost, as described in paragraph (b) above, is included with the
distribution. This option is valid only if you are a non-commercial
party, and only if you received the Object Form of the Work along with
such an offer.

You may copy and distribute the Work either gratis or for a fee, and
if desired, you may offer warranty protection for the Work.

The aggregation of the Work with other works that are not based on the
Work -- such as but not limited to inclusion in a publication,
broadcast, compilation, or other media -- does not bring the other
works in the scope of the License; nor does such aggregation void the
terms of the License for the Work.

4. MODIFICATION.

Permission is granted to modify or sample from a copy of the Work,
producing a derivative work, and to distribute the derivative work
under the terms described in the section for distribution above,
provided that the following terms are met:

(a) The new, derivative work is published under the terms of this
License.

(b) The derivative work is given a new name, so that its name or title
cannot be confused with the Work, or with a version of the Work, in
any way.

(c) Appropriate authorship credit is given: for the differences
between the Work and the new derivative work, authorship is attributed
to you, while the material sampled or used from the Work remains
attributed to the original Author; appropriate notice must be included
with the new work indicating the nature and the dates of any
modifications of the Work made by you.

5. NO RESTRICTIONS.

You may not impose any further restrictions on the Work or any of its
derivative works beyond those restrictions described in this License.

6. ACCEPTANCE.

Copying, distributing or modifying the Work (including but not limited
to sampling from the Work in a new work) indicates acceptance of these
terms. If you do not follow the terms of this License, any rights
granted to you by the License are null and void. The copying,
distribution or modification of the Work outside of the terms
described in this License is expressly prohibited by law.

If for any reason, conditions are imposed on you that forbid you to
fulfill the conditions of this License, you may not copy, distribute
or modify the Work at all.

If any part of this License is found to be in conflict with the law,
that part shall be interpreted in its broadest meaning consistent with
the law, and no other parts of the License shall be affected.

7. NO WARRANTY.

THE WORK IS PROVIDED "AS IS," AND COMES WITH ABSOLUTELY NO WARRANTY,
EXPRESS OR IMPLIED, TO THE EXTENT PERMITTED BY APPLICABLE LAW,
INCLUDING BUT NOT LIMITED TO THE IMPLIED WARRANTIES OF MERCHANTABILITY
OR FITNESS FOR A PARTICULAR PURPOSE.

8. DISCLAIMER OF LIABILITY.

IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
IN ANY WAY OUT OF THE USE OF THIS WORK, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.

END OF TERMS AND CONDITIONS


## Private decoy and unique-ID random helper: Boost Software License 1.0

The native MT19937-64 and bounded-integer mapping in
`src/chemistry/decoy_random.rs` are derived from the inspected Boost 1.90
`random/mersenne_twister.hpp` and `random/uniform_int_distribution.hpp`.
Unique ID generation reuses its raw engine without changing the recurrence.
The derived helper retains their component license and author notices:

Copyright Jens Maurer 2000-2001

Copyright Steven Watanabe 2010, 2011

Boost Software License - Version 1.0 - August 17th, 2003

Permission is hereby granted, free of charge, to any person or organization
obtaining a copy of the software and accompanying documentation covered by
this license (the "Software") to use, reproduce, display, distribute,
execute, and transmit the Software, and to prepare derivative works of the
Software, and to permit third-parties to whom the Software is furnished to
do so, all subject to the following:

The copyright notices in the Software and this entire statement, including
the above license grant, this restriction and the following disclaimer,
must be included in all copies of the Software, in whole or in part, and
all derivative works of the Software, unless such copies or derivative
works are solely in the form of machine-executable object code generated by
a source language processor.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.

The unchanged license text is available from the
[Boost license page](https://www.boost.org/LICENSE_1_0.txt).

## Native runtime clock dependencies

Local log timestamps use `chrono 0.4.45` with its clock feature. Process CPU
timing on Unix/Windows uses `cpu-time 1.0.0`. Both are licensed MIT OR Apache-2.0;
platform adapter dependencies and their exact versions are recorded in Cargo.lock.
Their source archives retain the individual dependency license notices.


## Raw Numpress codec: Johan Teleman BSD-3-Clause

The native port in `src/format/numpress.rs` derives from the pinned MSNumpress
implementation. Its original component notice is reproduced here for source
and binary redistribution. This is not a new runtime dependency.

```text
Native Rust MSNumpress codecs
        johan.teleman@immun.lth.se

        This distribution goes under the BSD 3-clause license. If you prefer to use Apache
        version 2.0, that is also available at https://github.com/fickludd/ms-numpress
        Copyright (c) 2013, Johan Teleman
        All rights reserved.

        Redistribution and use in source and binary forms, with or without modification,
        are permitted provided that the following conditions are met:

*         Redistributions of source code must retain the above copyright notice, this list
        of conditions and the following disclaimer.
*        Redistributions in binary form must reproduce the above copyright notice, this
        list of conditions and the following disclaimer in the documentation and/or other
        materials provided with the distribution.
*        Neither the name of the Lund University nor the names of its contributors may be
        used to endorse or promote products derived from this software without specific
        prior written permission.

        THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
        EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
        OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
        SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
        SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
        OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
        HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
        OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
        SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```
