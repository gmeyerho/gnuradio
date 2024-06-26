/* -*- c++ -*- */
/*
 * Copyright 2015 Free Software Foundation, Inc.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "ldpc_gen_mtrx_encoder_impl.h"

namespace gr {
namespace fec {
namespace code {

generic_encoder::sptr ldpc_gen_mtrx_encoder::make(const ldpc_G_matrix::sptr G_obj)
{
    return generic_encoder::sptr(new ldpc_gen_mtrx_encoder_impl(G_obj));
}

ldpc_gen_mtrx_encoder_impl::ldpc_gen_mtrx_encoder_impl(const ldpc_G_matrix::sptr G_obj)
    : generic_encoder("ldpc_gen_mtrx_encoder")
{
    // Generator matrix to use for encoding
    d_G = G_obj;

    d_rate = static_cast<double>(d_G->n()) / static_cast<double>(d_G->k());

    // Set frame size to k, the # of bits in the information word
    // All buffers and settings will be based on this value.
    set_frame_size(d_G->k());
}

ldpc_gen_mtrx_encoder_impl::~ldpc_gen_mtrx_encoder_impl() {}

int ldpc_gen_mtrx_encoder_impl::get_output_size() { return d_output_size; }

int ldpc_gen_mtrx_encoder_impl::get_input_size() { return d_frame_size; }

bool ldpc_gen_mtrx_encoder_impl::set_frame_size(unsigned int frame_size)
{
    bool ret = true;

    if (frame_size % d_G->k() != 0) {
        d_logger->error("Frame size ({:d} bits) must be a "
                        "multiple of the information word "
                        "size of the LDPC matrix ({:d}).",
                        frame_size,
                        d_G->k());
        throw std::runtime_error("ldpc_gen_mtrx_encoder: cannot use frame size.");
    }

    d_frame_size = frame_size;

    d_output_size = static_cast<int>(round(d_rate * d_frame_size));

    return ret;
}

double ldpc_gen_mtrx_encoder_impl::rate() { return d_rate; }

void ldpc_gen_mtrx_encoder_impl::generic_work(const void* inbuffer, void* outbuffer)
{
    // Populate the information word
    const unsigned char* in = (const unsigned char*)inbuffer;
    unsigned char* out = (unsigned char*)outbuffer;

    int j = 0;
    for (int i = 0; i < get_input_size(); i += d_G->k()) {
        d_G->encode(&out[j], &in[i]);
        j += d_G->n();
    }
}

} /* namespace code */
} /* namespace fec */
} /* namespace gr */
