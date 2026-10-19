# Library: https://github.com/tyuyu-62/cp-library
from array import array
from gc import collect


_NTT_998 = (
    911660635,
    86583718,
    (
        0, 911660635, 509520358, 369330050, 332049552, 983190778, 123842337, 238493703, 975955924,
        603855026, 856644456, 131300601, 842657263, 730768835, 942482514, 806263778, 151565301, 510815449,
        503497456, 743006876, 741047443, 56250497, 867605899, 0,
    ),
    (
        0, 86583718, 372528824, 373294451, 645684063, 112220581, 692852209, 155456985, 797128860, 90816748,
        860285882, 927414960, 354738543, 109331171, 293255632, 535113200, 308540755, 121186627, 608385704,
        438932459, 359477183, 824071951, 103369235, 0,
    ),
    (
        0, 372528824, 337190230, 454590761, 816400692, 578227951, 180142363, 83780245, 6597683, 70046822,
        623238099, 183021267, 402682409, 631680428, 344509872, 689220186, 365017329, 774342554, 729444058,
        102986190, 128751033, 395565204, 0, 0,
    ),
    (
        0, 509520358, 929031873, 170256584, 839780419, 282974284, 395914482, 444904435, 72135471, 638914820,
        66769500, 771127074, 985925487, 262319669, 262341272, 625870173, 768022760, 859816005, 914661783,
        430819711, 272774365, 530924681, 0, 0,
    ),
)
_MAX_NTT_SIZE = 1 << 23
_INVERSE_NTT_SIZE = 998244234

def _ntt998(a: list[int], invert: bool, half_input: bool = False) -> None:
    n = len(a)
    h = (n - 1).bit_length()
    (
        imag,
        inverse_imag,
        rate2,
        inverse_rate2,
        rate3,
        inverse_rate3,
    ) = _NTT_998
    if not invert:
        length = 0
        if half_input and n > 1:
            for i in range(n // 2):
                a[n // 2 + i] = a[i]
            length = 1
        while length < h:
            if h - length == 1:
                half = 1 << (h - length - 1)
                rotation = 1
                for block in range(1 << length):
                    offset = block << (h - length)
                    for i in range(half):
                        left = a[offset + i]
                        right = a[offset + half + i] * rotation
                        a[offset + i] = (left + right) % 998244353
                        a[offset + half + i] = (left - right) % 998244353
                    rotation = (
                        rotation
                        * rate2[(~block & -~block).bit_length()]
                        % 998244353
                    )
                length += 1
            else:
                quarter = 1 << (h - length - 2)
                rotation = 1
                for block in range(1 << length):
                    rotation2 = rotation * rotation % 998244353
                    rotation3 = rotation2 * rotation % 998244353
                    offset = block << (h - length)
                    for i in range(quarter):
                        a0 = a[offset + i]
                        a1 = a[offset + quarter + i] * rotation
                        a2 = a[offset + quarter * 2 + i] * rotation2
                        a3 = a[offset + quarter * 3 + i] * rotation3
                        imaginary = (a1 - a3) % 998244353 * imag
                        a[offset + i] = (a0 + a2 + a1 + a3) % 998244353
                        a[offset + quarter + i] = (
                            a0 + a2 - a1 - a3
                        ) % 998244353
                        a[offset + quarter * 2 + i] = (
                            a0 - a2 + imaginary
                        ) % 998244353
                        a[offset + quarter * 3 + i] = (
                            a0 - a2 - imaginary
                        ) % 998244353
                    rotation = (
                        rotation
                        * rate3[(~block & -~block).bit_length()]
                        % 998244353
                    )
                length += 2
        return

    length = h
    while length:
        if length == 1:
            half = 1 << (h - length)
            rotation = 1
            for block in range(1 << (length - 1)):
                offset = block << (h - length + 1)
                for i in range(half):
                    left = a[offset + i]
                    right = a[offset + half + i]
                    a[offset + i] = (left + right) % 998244353
                    a[offset + half + i] = (
                        (left - right) * rotation % 998244353
                    )
                rotation = (
                    rotation
                    * inverse_rate2[
                        (~block & -~block).bit_length()
                    ]
                    % 998244353
                )
            length -= 1
        else:
            quarter = 1 << (h - length)
            rotation = 1
            for block in range(1 << (length - 2)):
                rotation2 = rotation * rotation % 998244353
                rotation3 = rotation2 * rotation % 998244353
                offset = block << (h - length + 2)
                for i in range(quarter):
                    a0 = a[offset + i]
                    a1 = a[offset + quarter + i]
                    a2 = a[offset + quarter * 2 + i]
                    a3 = a[offset + quarter * 3 + i]
                    imaginary = (a2 - a3) * inverse_imag % 998244353
                    a[offset + i] = (a0 + a1 + a2 + a3) % 998244353
                    a[offset + quarter + i] = (
                        (a0 - a1 + imaginary) * rotation % 998244353
                    )
                    a[offset + quarter * 2 + i] = (
                        (a0 + a1 - a2 - a3) * rotation2 % 998244353
                    )
                    a[offset + quarter * 3 + i] = (
                        (a0 - a1 - imaginary) * rotation3 % 998244353
                    )
                rotation = (
                    rotation
                    * inverse_rate3[
                        (~block & -~block).bit_length()
                    ]
                    % 998244353
                )
            length -= 2


def _quarter_product(x, y, branch):
    h = (len(x) - 1).bit_length()
    levels = h - 3
    imag, inverse_imag, _, _, rate3, inverse_rate3 = _NTT_998
    root_power = (0, 2, 1, 3)[branch]
    forward_rotations = [pow(3, (998244352 >> (length + 4)) * root_power, 998244353) for length in range(0, levels, 2)]
    inverse_rotations = [pow(332748118, (998244352 >> (length + 2)) * root_power, 998244353) for length in range(2, levels + 1, 2)]
    weight = pow(3, (998244352 // (len(x) // 2)) * root_power, 998244353)
    rates = (0, 998244352, 86583718, 488723995, 628914303, 666194801, 15053575, 874402016, 759750650, 22288429, 394389327, 141599897, 866943752, 155587090, 267475518, 55761839, 191980575, 846679052, 487428904, 494746897, 255237477, 257196910, 941993856, 0)

    def forward(a, start, length, count):
        quarter = 1 << (h - length - 2)
        rotation = forward_rotations[length // 2]
        first = (branch << length) + (start >> (h - length))
        for block in range(count):
            rotation2 = rotation * rotation % 998244353
            rotation3 = rotation2 * rotation % 998244353
            offset = start + (block << (h - length))
            for i in range(quarter):
                a0 = a[offset + i]
                a1 = a[offset + quarter + i] * rotation
                a2 = a[offset + quarter * 2 + i] * rotation2
                a3 = a[offset + quarter * 3 + i] * rotation3
                imaginary = (a1 - a3) % 998244353 * imag
                a[offset + i] = (a0 + a2 + a1 + a3) % 998244353
                a[offset + quarter + i] = (
                    a0 + a2 - a1 - a3
                ) % 998244353
                a[offset + quarter * 2 + i] = (
                    a0 - a2 + imaginary
                ) % 998244353
                a[offset + quarter * 3 + i] = (
                    a0 - a2 - imaginary
                ) % 998244353
            rotation = rotation * rate3[(~(first + block) & -~(first + block)).bit_length()] % 998244353
        return rotation

    def inverse(start, length, count):
        quarter = 1 << (h - length)
        rotation = inverse_rotations[length // 2 - 1]
        first = (branch << (length - 2)) + (start >> (h - length + 2))
        a = x
        for block in range(count):
            rotation2 = rotation * rotation % 998244353
            rotation3 = rotation2 * rotation % 998244353
            offset = start + (block << (h - length + 2))
            for i in range(quarter):
                a0 = a[offset + i]
                a1 = a[offset + quarter + i]
                a2 = a[offset + quarter * 2 + i]
                a3 = a[offset + quarter * 3 + i]
                imaginary = (a2 - a3) * inverse_imag % 998244353
                a[offset + i] = (a0 + a1 + a2 + a3) % 998244353
                a[offset + quarter + i] = (
                    (a0 - a1 + imaginary) * rotation % 998244353
                )
                a[offset + quarter * 2 + i] = (
                    (a0 + a1 - a2 - a3) * rotation2 % 998244353
                )
                a[offset + quarter * 3 + i] = (
                    (a0 - a1 - imaginary) * rotation3 % 998244353
                )
            rotation = rotation * inverse_rate3[(~(first + block) & -~(first + block)).bit_length()] % 998244353
        inverse_rotations[length // 2 - 1] = rotation

    def bottom(start, count):
        length = levels - 2
        rotation = forward_rotations[length // 2]
        irotation = inverse_rotations[levels // 2 - 1]
        first = (branch << length) + (start >> 5)
        for block in range(count):
            offset = start + block * 32
            rotation2 = rotation * rotation % 998244353
            rotation3 = rotation2 * rotation % 998244353
            irotation2 = irotation * irotation % 998244353
            irotation3 = irotation2 * irotation % 998244353
            t0 = x[offset + 0]
            t1 = x[offset + 8] * rotation
            t2 = x[offset + 16] * rotation2
            t3 = x[offset + 24] * rotation3
            imaginary = (t1 - t3) % 998244353 * imag
            a0_0 = (t0 + t2 + t1 + t3) % 998244353
            a1_0 = (t0 + t2 - t1 - t3) % 998244353
            a2_0 = (t0 - t2 + imaginary) % 998244353
            a3_0 = (t0 - t2 - imaginary) % 998244353
            t0 = x[offset + 1]
            t1 = x[offset + 9] * rotation
            t2 = x[offset + 17] * rotation2
            t3 = x[offset + 25] * rotation3
            imaginary = (t1 - t3) % 998244353 * imag
            a0_1 = (t0 + t2 + t1 + t3) % 998244353
            a1_1 = (t0 + t2 - t1 - t3) % 998244353
            a2_1 = (t0 - t2 + imaginary) % 998244353
            a3_1 = (t0 - t2 - imaginary) % 998244353
            t0 = x[offset + 2]
            t1 = x[offset + 10] * rotation
            t2 = x[offset + 18] * rotation2
            t3 = x[offset + 26] * rotation3
            imaginary = (t1 - t3) % 998244353 * imag
            a0_2 = (t0 + t2 + t1 + t3) % 998244353
            a1_2 = (t0 + t2 - t1 - t3) % 998244353
            a2_2 = (t0 - t2 + imaginary) % 998244353
            a3_2 = (t0 - t2 - imaginary) % 998244353
            t0 = x[offset + 3]
            t1 = x[offset + 11] * rotation
            t2 = x[offset + 19] * rotation2
            t3 = x[offset + 27] * rotation3
            imaginary = (t1 - t3) % 998244353 * imag
            a0_3 = (t0 + t2 + t1 + t3) % 998244353
            a1_3 = (t0 + t2 - t1 - t3) % 998244353
            a2_3 = (t0 - t2 + imaginary) % 998244353
            a3_3 = (t0 - t2 - imaginary) % 998244353
            t0 = x[offset + 4]
            t1 = x[offset + 12] * rotation
            t2 = x[offset + 20] * rotation2
            t3 = x[offset + 28] * rotation3
            imaginary = (t1 - t3) % 998244353 * imag
            a0_4 = (t0 + t2 + t1 + t3) % 998244353
            a1_4 = (t0 + t2 - t1 - t3) % 998244353
            a2_4 = (t0 - t2 + imaginary) % 998244353
            a3_4 = (t0 - t2 - imaginary) % 998244353
            t0 = x[offset + 5]
            t1 = x[offset + 13] * rotation
            t2 = x[offset + 21] * rotation2
            t3 = x[offset + 29] * rotation3
            imaginary = (t1 - t3) % 998244353 * imag
            a0_5 = (t0 + t2 + t1 + t3) % 998244353
            a1_5 = (t0 + t2 - t1 - t3) % 998244353
            a2_5 = (t0 - t2 + imaginary) % 998244353
            a3_5 = (t0 - t2 - imaginary) % 998244353
            t0 = x[offset + 6]
            t1 = x[offset + 14] * rotation
            t2 = x[offset + 22] * rotation2
            t3 = x[offset + 30] * rotation3
            imaginary = (t1 - t3) % 998244353 * imag
            a0_6 = (t0 + t2 + t1 + t3) % 998244353
            a1_6 = (t0 + t2 - t1 - t3) % 998244353
            a2_6 = (t0 - t2 + imaginary) % 998244353
            a3_6 = (t0 - t2 - imaginary) % 998244353
            t0 = x[offset + 7]
            t1 = x[offset + 15] * rotation
            t2 = x[offset + 23] * rotation2
            t3 = x[offset + 31] * rotation3
            imaginary = (t1 - t3) % 998244353 * imag
            a0_7 = (t0 + t2 + t1 + t3) % 998244353
            a1_7 = (t0 + t2 - t1 - t3) % 998244353
            a2_7 = (t0 - t2 + imaginary) % 998244353
            a3_7 = (t0 - t2 - imaginary) % 998244353
            weight2 = rotation * imag % 998244353
            c0_0 = (a0_0 * a0_0 + rotation * ((2 * a0_1 * a0_7 + 2 * a0_2 * a0_6 + 2 * a0_3 * a0_5 + a0_4 * a0_4) % 998244353)) % 998244353
            c0_1 = (2 * a0_0 * a0_1 + rotation * ((2 * a0_2 * a0_7 + 2 * a0_3 * a0_6 + 2 * a0_4 * a0_5) % 998244353)) % 998244353
            c0_2 = (2 * a0_0 * a0_2 + a0_1 * a0_1 + rotation * ((2 * a0_3 * a0_7 + 2 * a0_4 * a0_6 + a0_5 * a0_5) % 998244353)) % 998244353
            c0_3 = (2 * a0_0 * a0_3 + 2 * a0_1 * a0_2 + rotation * ((2 * a0_4 * a0_7 + 2 * a0_5 * a0_6) % 998244353)) % 998244353
            c0_4 = (2 * a0_0 * a0_4 + 2 * a0_1 * a0_3 + a0_2 * a0_2 + rotation * ((2 * a0_5 * a0_7 + a0_6 * a0_6) % 998244353)) % 998244353
            c0_5 = (2 * a0_0 * a0_5 + 2 * a0_1 * a0_4 + 2 * a0_2 * a0_3 + rotation * ((2 * a0_6 * a0_7) % 998244353)) % 998244353
            c0_6 = (2 * a0_0 * a0_6 + 2 * a0_1 * a0_5 + 2 * a0_2 * a0_4 + a0_3 * a0_3 + rotation * ((a0_7 * a0_7) % 998244353)) % 998244353
            c0_7 = (2 * a0_0 * a0_7 + 2 * a0_1 * a0_6 + 2 * a0_2 * a0_5 + 2 * a0_3 * a0_4) % 998244353
            c1_0 = (a1_0 * a1_0 + (998244353 - rotation) * ((2 * a1_1 * a1_7 + 2 * a1_2 * a1_6 + 2 * a1_3 * a1_5 + a1_4 * a1_4) % 998244353)) % 998244353
            c1_1 = (2 * a1_0 * a1_1 + (998244353 - rotation) * ((2 * a1_2 * a1_7 + 2 * a1_3 * a1_6 + 2 * a1_4 * a1_5) % 998244353)) % 998244353
            c1_2 = (2 * a1_0 * a1_2 + a1_1 * a1_1 + (998244353 - rotation) * ((2 * a1_3 * a1_7 + 2 * a1_4 * a1_6 + a1_5 * a1_5) % 998244353)) % 998244353
            c1_3 = (2 * a1_0 * a1_3 + 2 * a1_1 * a1_2 + (998244353 - rotation) * ((2 * a1_4 * a1_7 + 2 * a1_5 * a1_6) % 998244353)) % 998244353
            c1_4 = (2 * a1_0 * a1_4 + 2 * a1_1 * a1_3 + a1_2 * a1_2 + (998244353 - rotation) * ((2 * a1_5 * a1_7 + a1_6 * a1_6) % 998244353)) % 998244353
            c1_5 = (2 * a1_0 * a1_5 + 2 * a1_1 * a1_4 + 2 * a1_2 * a1_3 + (998244353 - rotation) * ((2 * a1_6 * a1_7) % 998244353)) % 998244353
            c1_6 = (2 * a1_0 * a1_6 + 2 * a1_1 * a1_5 + 2 * a1_2 * a1_4 + a1_3 * a1_3 + (998244353 - rotation) * ((a1_7 * a1_7) % 998244353)) % 998244353
            c1_7 = (2 * a1_0 * a1_7 + 2 * a1_1 * a1_6 + 2 * a1_2 * a1_5 + 2 * a1_3 * a1_4) % 998244353
            c2_0 = (a2_0 * a2_0 + weight2 * ((2 * a2_1 * a2_7 + 2 * a2_2 * a2_6 + 2 * a2_3 * a2_5 + a2_4 * a2_4) % 998244353)) % 998244353
            c2_1 = (2 * a2_0 * a2_1 + weight2 * ((2 * a2_2 * a2_7 + 2 * a2_3 * a2_6 + 2 * a2_4 * a2_5) % 998244353)) % 998244353
            c2_2 = (2 * a2_0 * a2_2 + a2_1 * a2_1 + weight2 * ((2 * a2_3 * a2_7 + 2 * a2_4 * a2_6 + a2_5 * a2_5) % 998244353)) % 998244353
            c2_3 = (2 * a2_0 * a2_3 + 2 * a2_1 * a2_2 + weight2 * ((2 * a2_4 * a2_7 + 2 * a2_5 * a2_6) % 998244353)) % 998244353
            c2_4 = (2 * a2_0 * a2_4 + 2 * a2_1 * a2_3 + a2_2 * a2_2 + weight2 * ((2 * a2_5 * a2_7 + a2_6 * a2_6) % 998244353)) % 998244353
            c2_5 = (2 * a2_0 * a2_5 + 2 * a2_1 * a2_4 + 2 * a2_2 * a2_3 + weight2 * ((2 * a2_6 * a2_7) % 998244353)) % 998244353
            c2_6 = (2 * a2_0 * a2_6 + 2 * a2_1 * a2_5 + 2 * a2_2 * a2_4 + a2_3 * a2_3 + weight2 * ((a2_7 * a2_7) % 998244353)) % 998244353
            c2_7 = (2 * a2_0 * a2_7 + 2 * a2_1 * a2_6 + 2 * a2_2 * a2_5 + 2 * a2_3 * a2_4) % 998244353
            c3_0 = (a3_0 * a3_0 + (998244353 - weight2) * ((2 * a3_1 * a3_7 + 2 * a3_2 * a3_6 + 2 * a3_3 * a3_5 + a3_4 * a3_4) % 998244353)) % 998244353
            c3_1 = (2 * a3_0 * a3_1 + (998244353 - weight2) * ((2 * a3_2 * a3_7 + 2 * a3_3 * a3_6 + 2 * a3_4 * a3_5) % 998244353)) % 998244353
            c3_2 = (2 * a3_0 * a3_2 + a3_1 * a3_1 + (998244353 - weight2) * ((2 * a3_3 * a3_7 + 2 * a3_4 * a3_6 + a3_5 * a3_5) % 998244353)) % 998244353
            c3_3 = (2 * a3_0 * a3_3 + 2 * a3_1 * a3_2 + (998244353 - weight2) * ((2 * a3_4 * a3_7 + 2 * a3_5 * a3_6) % 998244353)) % 998244353
            c3_4 = (2 * a3_0 * a3_4 + 2 * a3_1 * a3_3 + a3_2 * a3_2 + (998244353 - weight2) * ((2 * a3_5 * a3_7 + a3_6 * a3_6) % 998244353)) % 998244353
            c3_5 = (2 * a3_0 * a3_5 + 2 * a3_1 * a3_4 + 2 * a3_2 * a3_3 + (998244353 - weight2) * ((2 * a3_6 * a3_7) % 998244353)) % 998244353
            c3_6 = (2 * a3_0 * a3_6 + 2 * a3_1 * a3_5 + 2 * a3_2 * a3_4 + a3_3 * a3_3 + (998244353 - weight2) * ((a3_7 * a3_7) % 998244353)) % 998244353
            c3_7 = (2 * a3_0 * a3_7 + 2 * a3_1 * a3_6 + 2 * a3_2 * a3_5 + 2 * a3_3 * a3_4) % 998244353
            imaginary = (c2_0 - c3_0) * inverse_imag % 998244353
            x[offset + 0] = (c0_0 + c1_0 + c2_0 + c3_0) % 998244353
            x[offset + 8] = ((c0_0 - c1_0 + imaginary) * irotation) % 998244353
            x[offset + 16] = ((c0_0 + c1_0 - c2_0 - c3_0) * irotation2) % 998244353
            x[offset + 24] = ((c0_0 - c1_0 - imaginary) * irotation3) % 998244353
            imaginary = (c2_1 - c3_1) * inverse_imag % 998244353
            x[offset + 1] = (c0_1 + c1_1 + c2_1 + c3_1) % 998244353
            x[offset + 9] = ((c0_1 - c1_1 + imaginary) * irotation) % 998244353
            x[offset + 17] = ((c0_1 + c1_1 - c2_1 - c3_1) * irotation2) % 998244353
            x[offset + 25] = ((c0_1 - c1_1 - imaginary) * irotation3) % 998244353
            imaginary = (c2_2 - c3_2) * inverse_imag % 998244353
            x[offset + 2] = (c0_2 + c1_2 + c2_2 + c3_2) % 998244353
            x[offset + 10] = ((c0_2 - c1_2 + imaginary) * irotation) % 998244353
            x[offset + 18] = ((c0_2 + c1_2 - c2_2 - c3_2) * irotation2) % 998244353
            x[offset + 26] = ((c0_2 - c1_2 - imaginary) * irotation3) % 998244353
            imaginary = (c2_3 - c3_3) * inverse_imag % 998244353
            x[offset + 3] = (c0_3 + c1_3 + c2_3 + c3_3) % 998244353
            x[offset + 11] = ((c0_3 - c1_3 + imaginary) * irotation) % 998244353
            x[offset + 19] = ((c0_3 + c1_3 - c2_3 - c3_3) * irotation2) % 998244353
            x[offset + 27] = ((c0_3 - c1_3 - imaginary) * irotation3) % 998244353
            imaginary = (c2_4 - c3_4) * inverse_imag % 998244353
            x[offset + 4] = (c0_4 + c1_4 + c2_4 + c3_4) % 998244353
            x[offset + 12] = ((c0_4 - c1_4 + imaginary) * irotation) % 998244353
            x[offset + 20] = ((c0_4 + c1_4 - c2_4 - c3_4) * irotation2) % 998244353
            x[offset + 28] = ((c0_4 - c1_4 - imaginary) * irotation3) % 998244353
            imaginary = (c2_5 - c3_5) * inverse_imag % 998244353
            x[offset + 5] = (c0_5 + c1_5 + c2_5 + c3_5) % 998244353
            x[offset + 13] = ((c0_5 - c1_5 + imaginary) * irotation) % 998244353
            x[offset + 21] = ((c0_5 + c1_5 - c2_5 - c3_5) * irotation2) % 998244353
            x[offset + 29] = ((c0_5 - c1_5 - imaginary) * irotation3) % 998244353
            imaginary = (c2_6 - c3_6) * inverse_imag % 998244353
            x[offset + 6] = (c0_6 + c1_6 + c2_6 + c3_6) % 998244353
            x[offset + 14] = ((c0_6 - c1_6 + imaginary) * irotation) % 998244353
            x[offset + 22] = ((c0_6 + c1_6 - c2_6 - c3_6) * irotation2) % 998244353
            x[offset + 30] = ((c0_6 - c1_6 - imaginary) * irotation3) % 998244353
            imaginary = (c2_7 - c3_7) * inverse_imag % 998244353
            x[offset + 7] = (c0_7 + c1_7 + c2_7 + c3_7) % 998244353
            x[offset + 15] = ((c0_7 - c1_7 + imaginary) * irotation) % 998244353
            x[offset + 23] = ((c0_7 + c1_7 - c2_7 - c3_7) * irotation2) % 998244353
            x[offset + 31] = ((c0_7 - c1_7 - imaginary) * irotation3) % 998244353
            index = (~(first + block) & -~(first + block)).bit_length()
            rotation = rotation * rate3[index] % 998244353
            irotation = irotation * inverse_rate3[index] % 998244353
        forward_rotations[length // 2] = rotation
        inverse_rotations[levels // 2 - 1] = irotation

    def visit(start, length):
        nonlocal weight
        size = 1 << (h - length)
        if size > 2048:
            rotation = forward(x, start, length, 1)
            if x is not y:
                forward(y, start, length, 1)
            forward_rotations[length // 2] = rotation
            quarter = size // 4
            for child in range(4):
                visit(start + child * quarter, length + 2)
            inverse(start, length + 2, 1)
            return
        if x is y:
            for stage in range(length, levels - 2, 2):
                forward_rotations[stage // 2] = forward(x, start, stage, 1 << (stage - length))
            bottom(start, size // 32)
            for stage in range(levels - 2, length, -2):
                inverse(start, stage, 1 << (stage - length - 2))
            return
        for stage in range(length, levels, 2):
            count = 1 << (stage - length)
            rotation = forward(x, start, stage, count)
            forward(y, start, stage, count)
            forward_rotations[stage // 2] = rotation
        first = branch * (len(x) // 8) + start // 8
        for j in range(size // 8):
            i = start + j * 8
            a0 = x[i + 0]
            a1 = x[i + 1]
            a2 = x[i + 2]
            a3 = x[i + 3]
            a4 = x[i + 4]
            a5 = x[i + 5]
            a6 = x[i + 6]
            a7 = x[i + 7]
            b0 = y[i + 0]
            b1 = y[i + 1]
            b2 = y[i + 2]
            b3 = y[i + 3]
            b4 = y[i + 4]
            b5 = y[i + 5]
            b6 = y[i + 6]
            b7 = y[i + 7]
            x[i + 0] = (a0 * b0 + weight * ((a1 * b7 + a2 * b6 + a3 * b5 + a4 * b4 + a5 * b3 + a6 * b2 + a7 * b1) % 998244353)) % 998244353
            x[i + 1] = (a0 * b1 + a1 * b0 + weight * ((a2 * b7 + a3 * b6 + a4 * b5 + a5 * b4 + a6 * b3 + a7 * b2) % 998244353)) % 998244353
            x[i + 2] = (a0 * b2 + a1 * b1 + a2 * b0 + weight * ((a3 * b7 + a4 * b6 + a5 * b5 + a6 * b4 + a7 * b3) % 998244353)) % 998244353
            x[i + 3] = (a0 * b3 + a1 * b2 + a2 * b1 + a3 * b0 + weight * ((a4 * b7 + a5 * b6 + a6 * b5 + a7 * b4) % 998244353)) % 998244353
            x[i + 4] = (a0 * b4 + a1 * b3 + a2 * b2 + a3 * b1 + a4 * b0 + weight * ((a5 * b7 + a6 * b6 + a7 * b5) % 998244353)) % 998244353
            x[i + 5] = (a0 * b5 + a1 * b4 + a2 * b3 + a3 * b2 + a4 * b1 + a5 * b0 + weight * ((a6 * b7 + a7 * b6) % 998244353)) % 998244353
            x[i + 6] = (a0 * b6 + a1 * b5 + a2 * b4 + a3 * b3 + a4 * b2 + a5 * b1 + a6 * b0 + weight * ((a7 * b7) % 998244353)) % 998244353
            x[i + 7] = (a0 * b7 + a1 * b6 + a2 * b5 + a3 * b4 + a4 * b3 + a5 * b2 + a6 * b1 + a7 * b0) % 998244353
            k = first + j
            weight = weight * rates[(~k & -~k).bit_length()] % 998244353
        for stage in range(levels, length, -2):
            inverse(start, stage, 1 << (stage - length - 2))

    visit(0, 0)


def _convolution_quarters(a, b):
    size = _MAX_NTT_SIZE
    length = len(a) + len(b) - 1
    rates = (0, 998244352, 86583718, 488723995, 628914303, 666194801, 15053575, 874402016, 759750650, 22288429, 394389327, 141599897, 866943752, 155587090, 267475518, 55761839, 191980575, 846679052, 487428904, 494746897, 255237477, 257196910, 941993856, 0)
    images = []
    for branch in range(4):
        factor = (1, 998244352, 911660635, 86583718)[branch]

        def forward(values):
            x = array('I', [0]) * size
            count = len(values) - size
            if branch == 0:
                for i in range(count):
                    x[i] = (values[i] + values[size + i]) % 998244353
            elif branch == 1:
                for i in range(count):
                    x[i] = (values[i] - values[size + i]) % 998244353
            else:
                for i in range(count):
                    x[i] = (values[i] + factor * values[size + i]) % 998244353
            for i in range(count, size):
                x[i] = values[i] % 998244353
            return x

        x = forward(a)
        y = x if a is b else forward(b)
        _quarter_product(x, y, branch)
        del y
        images.append(x)
        collect()
    result = [0] * length
    a0, a1, a2, a3 = images
    inverse = 2 * _INVERSE_NTT_SIZE % 998244353
    inverse_imaginary = 86583718 * inverse % 998244353
    for i in range(size):
        u0 = a0[i] + a1[i]
        u1 = (a0[i] - a1[i]) * inverse
        u2 = a2[i] + a3[i]
        u3 = (a2[i] - a3[i]) * inverse_imaginary
        result[i] = (u0 + u2) * inverse % 998244353
        result[size + i] = (u1 + u3) % 998244353
        if 2 * size + i < length:
            result[2 * size + i] = (u0 - u2) * inverse % 998244353
        if 3 * size + i < length:
            result[3 * size + i] = (u1 - u3) % 998244353
    del images, a0, a1, a2, a3, x
    collect()
    return result


def _convolution_cosets(a, b):
    size = _MAX_NTT_SIZE
    length = len(a) + len(b) - 1
    count = (length + size - 1) // size
    root = 779057549
    points = [pow(root, i * size, 998244353) for i in range(count)]
    if count > 2:
        matrix = [
            [pow(x, j, 998244353) for j in range(count)]
            + [int(i == j) for j in range(count)]
            for i, x in enumerate(points)
        ]
        for i in range(count):
            inverse = pow(matrix[i][i], -1, 998244353)
            matrix[i] = [v * inverse % 998244353 for v in matrix[i]]
            for j in range(count):
                if i != j:
                    factor = matrix[j][i]
                    matrix[j] = [
                        (u - factor * v) % 998244353
                        for u, v in zip(matrix[j], matrix[i])
                    ]
        weights = [row[count:] for row in matrix]
    images = []
    inverse_size = _INVERSE_NTT_SIZE
    for k in range(count):
        if k:
            powers = [pow(root, k * j, 998244353) for j in range(7)]

        def forward(values):
            x = array('I', [0]) * size
            if k == 0:
                for j in range(min(size, len(values))):
                    x[j] = values[j] % 998244353
                for j in range(size, len(values)):
                    i = j & (size - 1)
                    x[i] = (x[i] + values[j]) % 998244353
            else:
                for j in range(min(size, len(values))):
                    x[j] = values[j] * powers[j % 7] % 998244353
                for j in range(size, len(values)):
                    i = j & (size - 1)
                    x[i] = (x[i] + values[j] * powers[j % 7]) % 998244353
            _ntt998(x, False, len(values) <= size // 2)
            return x
        x = forward(a)
        if a is b:
            for i in range(size):
                x[i] = x[i] * x[i] % 998244353
        else:
            y = forward(b)
            for i in range(size):
                x[i] = x[i] * y[i] % 998244353
            del y
        _ntt998(x, True)
        if k == 0:
            for i in range(size):
                x[i] = x[i] * inverse_size % 998244353
        else:
            factors = [
                inverse_size * pow(root, -k * j, 998244353) % 998244353
                for j in range(7)
            ]
            for i in range(size):
                x[i] = x[i] * factors[i % 7] % 998244353
        images.append(x)
        collect()
    result = [0] * length
    if count == 2:
        a0, a1 = images
        inverse = pow(points[1] - 1, -1, 998244353)
        for i in range(size):
            upper = (a1[i] - a0[i]) * inverse % 998244353
            result[i] = (a0[i] - upper) % 998244353
            if size + i < length:
                result[size + i] = upper
        del a0, a1
    elif count == 3:
        a0, a1, a2 = images
        w00, w01, w02 = weights[0]
        w10, w11, w12 = weights[1]
        w20, w21, w22 = weights[2]
        for i in range(size):
            v0, v1, v2 = a0[i], a1[i], a2[i]
            result[i] = (w00 * v0 + w01 * v1 + w02 * v2) % 998244353
            result[size + i] = (w10 * v0 + w11 * v1 + w12 * v2) % 998244353
            if 2 * size + i < length:
                result[2 * size + i] = (w20 * v0 + w21 * v1 + w22 * v2) % 998244353
        del a0, a1, a2
    else:
        a0, a1, a2, a3 = images
        w00, w01, w02, w03 = weights[0]
        w10, w11, w12, w13 = weights[1]
        w20, w21, w22, w23 = weights[2]
        w30, w31, w32, w33 = weights[3]
        for i in range(size):
            v0, v1, v2, v3 = a0[i], a1[i], a2[i], a3[i]
            result[i] = (w00 * v0 + w01 * v1 + w02 * v2 + w03 * v3) % 998244353
            result[size + i] = (w10 * v0 + w11 * v1 + w12 * v2 + w13 * v3) % 998244353
            result[2 * size + i] = (w20 * v0 + w21 * v1 + w22 * v2 + w23 * v3) % 998244353
            if 3 * size + i < length:
                result[3 * size + i] = (w30 * v0 + w31 * v1 + w32 * v2 + w33 * v3) % 998244353
        del a0, a1, a2, a3
    del images, x
    collect()
    return result


def convolution(a: list[int], b: list[int]) -> list[int]:
    """a and b are coefficient sequences in increasing degree order. Return a new list of
    length len(a)+len(b)-1 with c[k]=sum(a[i]*b[k-i]) over valid indices modulo
    998244353. If either input is empty, return [].

    Constraints: Integer lists or arrays of lengths N,M.
    Time: O(S log S) for S=next_power_of_two(N+M-1)<=2^23; O(N*M) if min(N,M)<=32;
        otherwise O((P+Q)*L log L+P*Q*L), B=2^22, L=2*B, P=ceil(N/B), Q=ceil(M/B).
    Space: O(N+M+2^23) for long outputs; O(S) otherwise.
    """
    if not a or not b:
        return []
    length = len(a) + len(b) - 1
    if min(len(a), len(b)) <= 32:
        x = [(value + 499122176) % 998244353 - 499122176 for value in a]
        y = x if a is b else [(value + 499122176) % 998244353 - 499122176 for value in b]
        result = [0] * length
        for i, value in enumerate(x):
            for j, other in enumerate(y):
                result[i + j] += value * other
        for i in range(length):
            result[i] %= 998244353
        return result
    ntt_size = _MAX_NTT_SIZE
    if length <= ntt_size:
        size = 1 << (length - 1).bit_length()
        x = array('I', [0]) * size
        for i in range(len(a)):
            x[i] = a[i] % 998244353
        _ntt998(x, False, len(a) <= size // 2)
        if a is b:
            for i in range(size):
                x[i] = x[i] * x[i] % 998244353
        else:
            y = array('I', [0]) * size
            for i in range(len(b)):
                y[i] = b[i] % 998244353
            _ntt998(y, False, len(b) <= size // 2)
            for i in range(size):
                x[i] = x[i] * y[i] % 998244353
            del y
        _ntt998(x, True)
        inverse = pow(size, 998244351, 998244353)
        result = [0] * length
        for i in range(length):
            result[i] = x[i] * inverse % 998244353
        return result
    if ntt_size < min(len(a), len(b)) and max(len(a), len(b)) <= ntt_size * 2:
        return _convolution_quarters(a, b)
    block = ntt_size // 2
    if length <= ntt_size * 4:
        count = (length + ntt_size - 1) // ntt_size
        p = (len(a) + block - 1) // block
        if a is b:
            use_cosets = 2 * count < 3 * p - 1
        else:
            q = (len(b) + block - 1) // block
            use_cosets = 3 * count < 2 * (p + q) - 1
        if use_cosets:
            return _convolution_cosets(a, b)

    frequency = array('I', [0]) * ntt_size

    def forward(values):
        spectra = []
        for start in range(0, len(values), block):
            count = min(block, len(values) - start)
            for i in range(count):
                frequency[i] = values[start + i] % 998244353
            for i in range(count, block):
                frequency[i] = 0
            _ntt998(frequency, False, True)
            spectra.append(array('I', frequency))
        return spectra
    left = forward(a)
    right = left if a is b else forward(b)
    result = [0] * length
    inverse = _INVERSE_NTT_SIZE
    for k in range(len(left) + len(right) - 1):
        first = max(0, k - len(right) + 1)
        for i in range(first, min(len(left), k + 1)):
            j = k - i
            if a is b and i > j:
                continue
            x, y = (left[i], right[j])
            factor = 2 if a is b and i != j else 1
            if i == first:
                for t in range(ntt_size):
                    frequency[t] = factor * x[t] * y[t] % 998244353
            else:
                for t in range(ntt_size):
                    frequency[t] = (frequency[t] + factor * x[t] * y[t]) % 998244353
        _ntt998(frequency, True)
        offset = k * block
        for t in range(min(ntt_size - 1, length - offset)):
            result[offset + t] = (result[offset + t] + frequency[t] * inverse) % 998244353
    del left, right, x, y, frequency
    collect()
    return result
