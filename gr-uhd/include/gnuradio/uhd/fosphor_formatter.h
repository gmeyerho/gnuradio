/* -*- c++ -*- */
/*
 * Copyright 2023 Ettus Research, A National Instruments Brand
 *
 * This file is part of GNU Radio
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_GR_UHD_FOSPHOR_FORMATTER_H
#define INCLUDED_GR_UHD_FOSPHOR_FORMATTER_H

#include <gnuradio/hier_block2.h>
#include <gnuradio/uhd/api.h>

namespace gr {
namespace uhd {

/*!
 * \brief Generic CPM modulator
 * \ingroup modulators_blk
 */
class GR_UHD_API fosphor_formatter : virtual public hier_block2
{
public:
    using sptr = std::shared_ptr<fosphor_formatter>;

    ~fosphor_formatter() override {}

    /*!
     * Make formatter block for fosphor_display
     *
     */
    static sptr make(int fft_size,
                     int num_bins,
                     int input_decim,
                     int waterfall_decim,
                     int histo_decim,
                     double scale,
                     double alpha,
                     double epsilon,
                     double trise,
                     double tdecay);
};

} /* namespace uhd */
} /* namespace gr */

#endif /* INCLUDED_GR_UHD_FOSPHOR_FORMATTER_H */
