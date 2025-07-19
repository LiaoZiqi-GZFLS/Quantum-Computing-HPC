    Matrix(int idx)
    {
        switch (idx)
        {
        case 0:
            data[0][0] = 1.0;
            data[0][1] = 1.0;
            data[1][0] = 1.0;
            data[1][1] = -1.0;
            powOFsqrt2_inv = 1;
            break;
        case 1:
            data[0][0] = 0.0;
            data[0][1] = 1.0;
            data[1][0] = 1.0;
            data[1][1] = 0.0;
            powOFsqrt2_inv = 0;
            break;
        case 2:
            data[0][0] = 0.0;
            data[0][1] = std::complex<double>(0, -1.0);
            data[1][0] = std::complex<double>(0, 1.0);
            data[1][1] = 0.0;
            powOFsqrt2_inv = 0;
            break;
        case 3:
            data[0][0] = 1.0;
            data[0][1] = 0.0;
            data[1][0] = 0.0;
            data[1][1] = -1.0;
            powOFsqrt2_inv = 0;
            break;
        case 4:
            data[0][0] = 1.0;
            data[0][1] = 0.0;
            data[1][0] = 0.0;
            data[1][1] = std::complex<double>(0, 1.0);
            powOFsqrt2_inv = 0;
            break;
        case 5:
            data[0][0] = std::complex<double>(0, 1.0);
            data[0][1] = std::complex<double>(0, -1.0);
            data[1][0] = std::complex<double>(0, 1.0);
            data[1][1] = std::complex<double>(0, 1.0);
            powOFsqrt2_inv = 1;
            break;
        case 6:
            data[0][0] = 1.0;
            data[0][1] = -1.0;
            data[1][0] = -1.0;
            data[1][1] = -1.0;
            powOFsqrt2_inv = 1;
            break;
        case 7:
            data[0][0] = 1.0;
            data[0][1] = std::complex<double>(0, 1.0);
            data[1][0] = -1.0;
            data[1][1] = std::complex<double>(0, 1.0);
            powOFsqrt2_inv = 1;
            break;
        case 8:
            data[0][0] = 0.0;
            data[0][1] = std::complex<double>(0, 1.0);
            data[1][0] = std::complex<double>(0, -1.0);
            data[1][1] = 0.0;
            powOFsqrt2_inv = 0;
            break;
        case 9:
            data[0][0] = std::complex<double>(0, -1.0);
            data[0][1] = std::complex<double>(0, 1.0);
            data[1][0] = std::complex<double>(0, -1.0);
            data[1][1] = std::complex<double>(0, -1.0);
            powOFsqrt2_inv = 1;
            break;
        case 10:
            data[0][0] = std::complex<double>(0, 1.0);
            data[0][1] = std::complex<double>(0, 1.0);
            data[1][0] = std::complex<double>(0, -1.0);
            data[1][1] = std::complex<double>(0, 1.0);
            powOFsqrt2_inv = 1;
            break;
        case 11:
            data[0][0] = std::complex<double>(0, 1.0);
            data[0][1] = 1.0;
            data[1][0] = std::complex<double>(0, -1.0);
            data[1][1] = 1.0;
            powOFsqrt2_inv = 1;
            break;
        case 12:
            data[0][0] = -1.0;
            data[0][1] = 1.0;
            data[1][0] = 1.0;
            data[1][1] = 1.0;
            powOFsqrt2_inv = 1;
            break;
        case 13:
            data[0][0] = std::complex<double>(0, -1.0);
            data[0][1] = std::complex<double>(0, -1.0);
            data[1][0] = std::complex<double>(0, 1.0);
            data[1][1] = std::complex<double>(0, -1.0);
            powOFsqrt2_inv = 1;
            break;
        case 14:
            data[0][0] = 1.0;
            data[0][1] = std::complex<double>(0, -1.0);
            data[1][0] = 1.0;
            data[1][1] = std::complex<double>(0, 1.0);
            powOFsqrt2_inv = 1;
            break;
        case 15:
            data[0][0] = std::complex<double>(1.0, 1.0);
            data[0][1] = std::complex<double>(1.0, -1.0);
            data[1][0] = std::complex<double>(1.0, -1.0);
            data[1][1] = std::complex<double>(1.0, 1.0);
            powOFsqrt2_inv = 0;
            break;
        case 16:
            data[0][0] = -1.0;
            data[0][1] = std::complex<double>(0, -1.0);
            data[1][0] = 1.0;
            data[1][1] = std::complex<double>(0, -1.0);
            powOFsqrt2_inv = 1;
            break;
        case 17:
            data[0][0] = 1.0;
            data[0][1] = -1.0;
            data[1][0] = 1.0;
            data[1][1] = 1.0;
            powOFsqrt2_inv = 1;
            break;
        case 18:
            data[0][0] = std::complex<double>(0, 1.0);
            data[0][1] = 0.0;
            data[1][0] = 0.0;
            data[1][1] = std::complex<double>(0, 1.0);
            powOFsqrt2_inv = 0;
            break;
        case 19:
            data[0][0] = std::complex<double>(0, 1.0);
            data[0][1] = 0.0;
            data[1][0] = 0.0;
            data[1][1] = 1.0;
            powOFsqrt2_inv = 0;
            break;
        case 20:
            data[0][0] = -1.0;
            data[0][1] = 0.0;
            data[1][0] = 0.0;
            data[1][1] = 1.0;
            powOFsqrt2_inv = 0;
            break;
        case 21:
            data[0][0] = std::complex<double>(0, -1.0);
            data[0][1] = 0.0;
            data[1][0] = 0.0;
            data[1][1] = std::complex<double>(0, -1.0);
            powOFsqrt2_inv = 0;
            break;
        case 22:
            data[0][0] = 0.0;
            data[0][1] = std::complex<double>(0, -1.0);
            data[1][0] = 1.0;
            data[1][1] = 0.0;
            powOFsqrt2_inv = 0;
            break;
        case 23:
            data[0][0] = std::complex<double>(0, 1.0);
            data[0][1] = std::complex<double>(0, -1.0);
            data[1][0] = 1.0;
            data[1][1] = 1.0;
            powOFsqrt2_inv = 1;
            break;
        case 24:
            data[0][0] = -1.0;
            data[0][1] = 0.0;
            data[1][0] = 0.0;
            data[1][1] = std::complex<double>(0, -1.0);
            powOFsqrt2_inv = 0;
            break;
        case 25:
            data[0][0] = 0.0;
            data[0][1] = -1.0;
            data[1][0] = 1.0;
            data[1][1] = 0.0;
            powOFsqrt2_inv = 0;
            break;
        case 26:
            data[0][0] = -1.0;
            data[0][1] = -1.0;
            data[1][0] = -1.0;
            data[1][1] = 1.0;
            powOFsqrt2_inv = 1;
            break;
        case 27:
            data[0][0] = std::complex<double>(0, -1.0);
            data[0][1] = -1.0;
            data[1][0] = std::complex<double>(0, 1.0);
            data[1][1] = -1.0;
            powOFsqrt2_inv = 1;
            break;
        case 28:
            data[0][0] = 0.0;
            data[0][1] = -1.0;
            data[1][0] = -1.0;
            data[1][1] = 0.0;
            powOFsqrt2_inv = 0;
            break;
        case 29:
            data[0][0] = std::complex<double>(0, -1.0);
            data[0][1] = 0.0;
            data[1][0] = 0.0;
            data[1][1] = -1.0;
            powOFsqrt2_inv = 0;
            break;
        case 30:
            data[0][0] = 0.0;
            data[0][1] = -1.0;
            data[1][0] = std::complex<double>(0, 1.0);
            data[1][1] = 0.0;
            powOFsqrt2_inv = 0;
            break;
        case 31:
            data[0][0] = 1.0;
            data[0][1] = -1.0;
            data[1][0] = std::complex<double>(0, 1.0);
            data[1][1] = std::complex<double>(0, 1.0);
            powOFsqrt2_inv = 1;
            break;
        case 32:
            data[0][0] = 0.0;
            data[0][1] = std::complex<double>(0, 1.0);
            data[1][0] = std::complex<double>(0, 1.0);
            data[1][1] = 0.0;
            powOFsqrt2_inv = 0;
            break;
        case 33:
            data[0][0] = 0.0;
            data[0][1] = std::complex<double>(0, 1.0);
            data[1][0] = -1.0;
            data[1][1] = 0.0;
            powOFsqrt2_inv = 0;
            break;
        case 34:
            data[0][0] = 0.0;
            data[0][1] = 1.0;
            data[1][0] = std::complex<double>(0, -1.0);
            data[1][1] = 0.0;
            powOFsqrt2_inv = 0;
            break;
        case 35:
            data[0][0] = 1.0;
            data[0][1] = 1.0;
            data[1][0] = std::complex<double>(0, -1.0);
            data[1][1] = std::complex<double>(0, 1.0);
            powOFsqrt2_inv = 1;
            break;
        case 36:
            data[0][0] = 1.0;
            data[0][1] = 0.0;
            data[1][0] = 0.0;
            data[1][1] = 1.0;
            powOFsqrt2_inv = 0;
            break;
        case 37:
            data[0][0] = 1.0;
            data[0][1] = std::complex<double>(0, 1.0);
            data[1][0] = std::complex<double>(0, 1.0);
            data[1][1] = 1.0;
            powOFsqrt2_inv = 1;
            break;
        case 38:
            data[0][0] = std::complex<double>(0, -1.0);
            data[0][1] = std::complex<double>(0, 1.0);
            data[1][0] = -1.0;
            data[1][1] = -1.0;
            powOFsqrt2_inv = 1;
            break;
        case 39:
            data[0][0] = 0.0;
            data[0][1] = 1.0;
            data[1][0] = -1.0;
            data[1][1] = 0.0;
            powOFsqrt2_inv = 0;
            break;
        case 40:
            data[0][0] = 1.0;
            data[0][1] = 1.0;
            data[1][0] = -1.0;
            data[1][1] = 1.0;
            powOFsqrt2_inv = 1;
            break;
        case 41:
            data[0][0] = 0.0;
            data[0][1] = std::complex<double>(0, -1.0);
            data[1][0] = std::complex<double>(0, -1.0);
            data[1][1] = 0.0;
            powOFsqrt2_inv = 0;
            break;
        case 42:
            data[0][0] = 1.0;
            data[0][1] = 0.0;
            data[1][0] = 0.0;
            data[1][1] = std::complex<double>(0, -1.0);
            powOFsqrt2_inv = 0;
            break;
        case 43:
            data[0][0] = std::complex<double>(0, 1.0);
            data[0][1] = -1.0;
            data[1][0] = std::complex<double>(0, 1.0);
            data[1][1] = 1.0;
            powOFsqrt2_inv = 1;
            break;
        case 44:
            data[0][0] = std::complex<double>(1.0, 1.0);
            data[0][1] = std::complex<double>(-1.0, 1.0);
            data[1][0] = std::complex<double>(-1.0, 1.0);
            data[1][1] = std::complex<double>(1.0, 1.0);
            powOFsqrt2_inv = 0;
            break;
        case 45:
            data[0][0] = std::complex<double>(0, 1.0);
            data[0][1] = std::complex<double>(0, 1.0);
            data[1][0] = std::complex<double>(0, 1.0);
            data[1][1] = std::complex<double>(0, -1.0);
            powOFsqrt2_inv = 1;
            break;
        case 46:
            data[0][0] = -1.0;
            data[0][1] = std::complex<double>(0, 1.0);
            data[1][0] = -1.0;
            data[1][1] = std::complex<double>(0, -1.0);
            powOFsqrt2_inv = 1;
            break;
        case 47:
            data[0][0] = std::complex<double>(0, -1.0);
            data[0][1] = 1.0;
            data[1][0] = std::complex<double>(0, -1.0);
            data[1][1] = -1.0;
            powOFsqrt2_inv = 1;
            break;
        case 48:
            data[0][0] = std::complex<double>(1.0, -1.0);
            data[0][1] = std::complex<double>(1.0, 1.0);
            data[1][0] = std::complex<double>(-1.0, -1.0);
            data[1][1] = std::complex<double>(-1.0, 1.0);
            powOFsqrt2_inv = 0;
            break;
        case 49:
            data[0][0] = std::complex<double>(1.0, 1.0);
            data[0][1] = std::complex<double>(1.0, 1.0);
            data[1][0] = std::complex<double>(-1.0, 1.0);
            data[1][1] = std::complex<double>(1.0, -1.0);
            powOFsqrt2_inv = 0;
            break;
        case 50:
            data[0][0] = std::complex<double>(-1.0, -1.0);
            data[0][1] = std::complex<double>(-1.0, 1.0);
            data[1][0] = std::complex<double>(-1.0, 1.0);
            data[1][1] = std::complex<double>(-1.0, -1.0);
            powOFsqrt2_inv = 0;
            break;
        case 51:
            data[0][0] = -1.0;
            data[0][1] = 1.0;
            data[1][0] = -1.0;
            data[1][1] = -1.0;
            powOFsqrt2_inv = 1;
            break;
        case 52:
            data[0][0] = std::complex<double>(0, -1.0);
            data[0][1] = std::complex<double>(0, -1.0);
            data[1][0] = std::complex<double>(0, -1.0);
            data[1][1] = std::complex<double>(0, 1.0);
            powOFsqrt2_inv = 1;
            break;
        case 53:
            data[0][0] = 1.0;
            data[0][1] = std::complex<double>(0, -1.0);
            data[1][0] = -1.0;
            data[1][1] = std::complex<double>(0, -1.0);
            powOFsqrt2_inv = 1;
            break;
        case 54:
            data[0][0] = -1.0;
            data[0][1] = 1.0;
            data[1][0] = std::complex<double>(0, -1.0);
            data[1][1] = std::complex<double>(0, -1.0);
            powOFsqrt2_inv = 1;
            break;
        case 55:
            data[0][0] = std::complex<double>(-1.0, -1.0);
            data[0][1] = std::complex<double>(1.0, -1.0);
            data[1][0] = std::complex<double>(1.0, -1.0);
            data[1][1] = std::complex<double>(-1.0, -1.0);
            powOFsqrt2_inv = 0;
            break;
        case 56:
            data[0][0] = std::complex<double>(-1.0, 1.0);
            data[0][1] = std::complex<double>(1.0, 1.0);
            data[1][0] = std::complex<double>(-1.0, -1.0);
            data[1][1] = std::complex<double>(1.0, -1.0);
            powOFsqrt2_inv = 0;
            break;
        case 57:
            data[0][0] = std::complex<double>(0, 1.0);
            data[0][1] = std::complex<double>(0, -1.0);
            data[1][0] = std::complex<double>(0, -1.0);
            data[1][1] = std::complex<double>(0, -1.0);
            powOFsqrt2_inv = 1;
            break;
        case 58:
            data[0][0] = std::complex<double>(1.0, 1.0);
            data[0][1] = std::complex<double>(-1.0, -1.0);
            data[1][0] = std::complex<double>(1.0, -1.0);
            data[1][1] = std::complex<double>(1.0, -1.0);
            powOFsqrt2_inv = 0;
            break;
        case 59:
            data[0][0] = std::complex<double>(0, 1.0);
            data[0][1] = 0.0;
            data[1][0] = 0.0;
            data[1][1] = std::complex<double>(0, -1.0);
            powOFsqrt2_inv = 0;
            break;
        case 60:
            data[0][0] = std::complex<double>(0, 1.0);
            data[0][1] = -1.0;
            data[1][0] = std::complex<double>(0, -1.0);
            data[1][1] = -1.0;
            powOFsqrt2_inv = 1;
            break;
        case 61:
            data[0][0] = std::complex<double>(-1.0, 1.0);
            data[0][1] = std::complex<double>(-1.0, -1.0);
            data[1][0] = std::complex<double>(1.0, 1.0);
            data[1][1] = std::complex<double>(1.0, -1.0);
            powOFsqrt2_inv = 0;
            break;
        case 62:
            data[0][0] = -1.0;
            data[0][1] = -1.0;
            data[1][0] = 1.0;
            data[1][1] = -1.0;
            powOFsqrt2_inv = 1;
            break;
        case 63:
            data[0][0] = std::complex<double>(1.0, -1.0);
            data[0][1] = std::complex<double>(-1.0, -1.0);
            data[1][0] = std::complex<double>(1.0, 1.0);
            data[1][1] = std::complex<double>(-1.0, 1.0);
            powOFsqrt2_inv = 0;
            break;
        case 64:
            data[0][0] = std::complex<double>(0, -1.0);
            data[0][1] = std::complex<double>(0, 1.0);
            data[1][0] = std::complex<double>(0, 1.0);
            data[1][1] = std::complex<double>(0, 1.0);
            powOFsqrt2_inv = 1;
            break;
        case 65:
            data[0][0] = std::complex<double>(1.0, -1.0);
            data[0][1] = std::complex<double>(-1.0, 1.0);
            data[1][0] = std::complex<double>(1.0, 1.0);
            data[1][1] = std::complex<double>(1.0, 1.0);
            powOFsqrt2_inv = 0;
            break;
        case 66:
            data[0][0] = 1.0;
            data[0][1] = std::complex<double>(0, 1.0);
            data[1][0] = 1.0;
            data[1][1] = std::complex<double>(0, -1.0);
            powOFsqrt2_inv = 1;
            break;
        case 67:
            data[0][0] = std::complex<double>(1.0, 1.0);
            data[0][1] = 0.0;
            data[1][0] = 0.0;
            data[1][1] = std::complex<double>(1.0, -1.0);
            powOFsqrt2_inv = 1;
            break;
        case 68:
            data[0][0] = std::complex<double>(-1.0, 1.0);
            data[0][1] = std::complex<double>(1.0, -1.0);
            data[1][0] = std::complex<double>(-1.0, -1.0);
            data[1][1] = std::complex<double>(-1.0, -1.0);
            powOFsqrt2_inv = 0;
            break;
        case 69:
            data[0][0] = std::complex<double>(1.0, 1.0);
            data[0][1] = std::complex<double>(-1.0, 1.0);
            data[1][0] = std::complex<double>(1.0, -1.0);
            data[1][1] = std::complex<double>(-1.0, -1.0);
            powOFsqrt2_inv = 0;
            break;
        case 70:
            data[0][0] = std::complex<double>(-1.0, -1.0);
            data[0][1] = std::complex<double>(-1.0, -1.0);
            data[1][0] = std::complex<double>(1.0, -1.0);
            data[1][1] = std::complex<double>(-1.0, 1.0);
            powOFsqrt2_inv = 0;
            break;
        case 71:
            data[0][0] = -1.0;
            data[0][1] = std::complex<double>(0, 1.0);
            data[1][0] = 1.0;
            data[1][1] = std::complex<double>(0, 1.0);
            powOFsqrt2_inv = 1;
            break;
        case 72:
            data[0][0] = 0.0;
            data[0][1] = std::complex<double>(0, 1.0);
            data[1][0] = 1.0;
            data[1][1] = 0.0;
            powOFsqrt2_inv = 0;
            break;
        case 73:
            data[0][0] = std::complex<double>(0, -1.0);
            data[0][1] = 0.0;
            data[1][0] = 0.0;
            data[1][1] = std::complex<double>(0, 1.0);
            powOFsqrt2_inv = 0;
            break;
        case 74:
            data[0][0] = std::complex<double>(0, -1.0);
            data[0][1] = 1.0;
            data[1][0] = std::complex<double>(0, 1.0);
            data[1][1] = 1.0;
            powOFsqrt2_inv = 1;
            break;
        case 75:
            data[0][0] = std::complex<double>(1.0, -1.0);
            data[0][1] = std::complex<double>(1.0, 1.0);
            data[1][0] = std::complex<double>(1.0, 1.0);
            data[1][1] = std::complex<double>(1.0, -1.0);
            powOFsqrt2_inv = 0;
            break;
        case 76:
            data[0][0] = std::complex<double>(0, 1.0);
            data[0][1] = std::complex<double>(0, 1.0);
            data[1][0] = -1.0;
            data[1][1] = 1.0;
            powOFsqrt2_inv = 1;
            break;
        case 77:
            data[0][0] = std::complex<double>(0, 1.0);
            data[0][1] = -1.0;
            data[1][0] = 1.0;
            data[1][1] = std::complex<double>(0, -1.0);
            powOFsqrt2_inv = 1;
            break;
        case 78:
            data[0][0] = std::complex<double>(0, 1.0);
            data[0][1] = 0.0;
            data[1][0] = 0.0;
            data[1][1] = -1.0;
            powOFsqrt2_inv = 0;
            break;
        case 79:
            data[0][0] = -1.0;
            data[0][1] = -1.0;
            data[1][0] = std::complex<double>(0, 1.0);
            data[1][1] = std::complex<double>(0, -1.0);
            powOFsqrt2_inv = 1;
            break;
        case 80:
            data[0][0] = -1.0;
            data[0][1] = 0.0;
            data[1][0] = 0.0;
            data[1][1] = -1.0;
            powOFsqrt2_inv = 0;
            break;
        case 81:
            data[0][0] = std::complex<double>(0, -1.0);
            data[0][1] = std::complex<double>(0, -1.0);
            data[1][0] = 1.0;
            data[1][1] = -1.0;
            powOFsqrt2_inv = 1;
            break;
        case 82:
            data[0][0] = std::complex<double>(0, -1.0);
            data[0][1] = -1.0;
            data[1][0] = 1.0;
            data[1][1] = std::complex<double>(0, 1.0);
            powOFsqrt2_inv = 1;
            break;
        case 83:
            data[0][0] = std::complex<double>(1.0, 1.0);
            data[0][1] = std::complex<double>(-1.0, 1.0);
            data[1][0] = std::complex<double>(1.0, 1.0);
            data[1][1] = std::complex<double>(1.0, -1.0);
            powOFsqrt2_inv = 0;
            break;
        case 84:
            data[0][0] = std::complex<double>(0, 1.0);
            data[0][1] = 1.0;
            data[1][0] = -1.0;
            data[1][1] = std::complex<double>(0, -1.0);
            powOFsqrt2_inv = 1;
            break;
        case 85:
            data[0][0] = std::complex<double>(0, 1.0);
            data[0][1] = std::complex<double>(0, 1.0);
            data[1][0] = 1.0;
            data[1][1] = -1.0;
            powOFsqrt2_inv = 1;
            break;
        case 86:
            data[0][0] = -1.0;
            data[0][1] = std::complex<double>(0, -1.0);
            data[1][0] = std::complex<double>(0, -1.0);
            data[1][1] = -1.0;
            powOFsqrt2_inv = 1;
            break;
        case 87:
            data[0][0] = -1.0;
            data[0][1] = 0.0;
            data[1][0] = 0.0;
            data[1][1] = std::complex<double>(0, 1.0);
            powOFsqrt2_inv = 0;
            break;
        case 88:
            data[0][0] = std::complex<double>(0, -1.0);
            data[0][1] = 0.0;
            data[1][0] = 0.0;
            data[1][1] = 1.0;
            powOFsqrt2_inv = 0;
            break;
        case 89:
            data[0][0] = std::complex<double>(0, -1.0);
            data[0][1] = std::complex<double>(0, 1.0);
            data[1][0] = 1.0;
            data[1][1] = 1.0;
            powOFsqrt2_inv = 1;
            break;
        case 90:
            data[0][0] = std::complex<double>(-1.0, -1.0);
            data[0][1] = std::complex<double>(1.0, 1.0);
            data[1][0] = std::complex<double>(-1.0, 1.0);
            data[1][1] = std::complex<double>(-1.0, 1.0);
            powOFsqrt2_inv = 0;
            break;
        case 91:
            data[0][0] = std::complex<double>(0, -1.0);
            data[0][1] = 1.0;
            data[1][0] = -1.0;
            data[1][1] = std::complex<double>(0, 1.0);
            powOFsqrt2_inv = 1;
            break;
        case 92:
            data[0][0] = -1.0;
            data[0][1] = std::complex<double>(0, 1.0);
            data[1][0] = std::complex<double>(0, 1.0);
            data[1][1] = -1.0;
            powOFsqrt2_inv = 1;
            break;
        case 93:
            data[0][0] = 0.0;
            data[0][1] = 1.0;
            data[1][0] = std::complex<double>(0, 1.0);
            data[1][1] = 0.0;
            powOFsqrt2_inv = 0;
            break;
        case 94:
            data[0][0] = std::complex<double>(1.0, -1.0);
            data[0][1] = std::complex<double>(1.0, 1.0);
            data[1][0] = std::complex<double>(-1.0, 1.0);
            data[1][1] = std::complex<double>(1.0, 1.0);
            powOFsqrt2_inv = 0;
            break;
        case 95:
            data[0][0] = 1.0;
            data[0][1] = std::complex<double>(0, -1.0);
            data[1][0] = std::complex<double>(0, -1.0);
            data[1][1] = 1.0;
            powOFsqrt2_inv = 1;
            break;
        case 96:
            data[0][0] = 1.0;
            data[0][1] = 1.0;
            data[1][0] = std::complex<double>(0, 1.0);
            data[1][1] = std::complex<double>(0, -1.0);
            powOFsqrt2_inv = 1;
            break;
        case 97:
            data[0][0] = std::complex<double>(0, 1.0);
            data[0][1] = 1.0;
            data[1][0] = std::complex<double>(0, 1.0);
            data[1][1] = -1.0;
            powOFsqrt2_inv = 1;
            break;
        case 98:
            data[0][0] = -1.0;
            data[0][1] = 1.0;
            data[1][0] = std::complex<double>(0, 1.0);
            data[1][1] = std::complex<double>(0, 1.0);
            powOFsqrt2_inv = 1;
            break;
        case 99:
            data[0][0] = 0.0;
            data[0][1] = std::complex<double>(0, -1.0);
            data[1][0] = -1.0;
            data[1][1] = 0.0;
            powOFsqrt2_inv = 0;
            break;
        case 100:
            data[0][0] = 0.0;
            data[0][1] = -1.0;
            data[1][0] = std::complex<double>(0, -1.0);
            data[1][1] = 0.0;
            powOFsqrt2_inv = 0;
            break;
        case 101:
            data[0][0] = std::complex<double>(1.0, 1.0);
            data[0][1] = std::complex<double>(1.0, -1.0);
            data[1][0] = std::complex<double>(-1.0, -1.0);
            data[1][1] = std::complex<double>(1.0, -1.0);
            powOFsqrt2_inv = 0;
            break;
        case 102:
            data[0][0] = 1.0;
            data[0][1] = -1.0;
            data[1][0] = std::complex<double>(0, -1.0);
            data[1][1] = std::complex<double>(0, -1.0);
            powOFsqrt2_inv = 1;
            break;
        case 103:
            data[0][0] = std::complex<double>(1.0, 1.0);
            data[0][1] = std::complex<double>(1.0, 1.0);
            data[1][0] = std::complex<double>(1.0, 1.0);
            data[1][1] = std::complex<double>(-1.0, -1.0);
            powOFsqrt2_inv = 0;
            break;
        case 104:
            data[0][0] = std::complex<double>(-1.0, -1.0);
            data[0][1] = std::complex<double>(-1.0, 1.0);
            data[1][0] = std::complex<double>(1.0, 1.0);
            data[1][1] = std::complex<double>(-1.0, 1.0);
            powOFsqrt2_inv = 0;
            break;
        case 105:
            data[0][0] = std::complex<double>(0, -1.0);
            data[0][1] = std::complex<double>(0, -1.0);
            data[1][0] = -1.0;
            data[1][1] = 1.0;
            powOFsqrt2_inv = 1;
            break;
        case 106:
            data[0][0] = 1.0;
            data[0][1] = std::complex<double>(0, -1.0);
            data[1][0] = std::complex<double>(0, 1.0);
            data[1][1] = -1.0;
            powOFsqrt2_inv = 1;
            break;
        case 107:
            data[0][0] = std::complex<double>(-1.0, -1.0);
            data[0][1] = std::complex<double>(1.0, -1.0);
            data[1][0] = std::complex<double>(-1.0, -1.0);
            data[1][1] = std::complex<double>(-1.0, 1.0);
            powOFsqrt2_inv = 0;
            break;
        case 108:
            data[0][0] = std::complex<double>(0, 1.0);
            data[0][1] = std::complex<double>(0, -1.0);
            data[1][0] = -1.0;
            data[1][1] = -1.0;
            powOFsqrt2_inv = 1;
            break;
        case 109:
            data[0][0] = std::complex<double>(1.0, 1.0);
            data[0][1] = std::complex<double>(1.0, -1.0);
            data[1][0] = std::complex<double>(-1.0, 1.0);
            data[1][1] = std::complex<double>(-1.0, -1.0);
            powOFsqrt2_inv = 0;
            break;
        case 110:
            data[0][0] = -1.0;
            data[0][1] = std::complex<double>(0, -1.0);
            data[1][0] = -1.0;
            data[1][1] = std::complex<double>(0, 1.0);
            powOFsqrt2_inv = 1;
            break;
        case 111:
            data[0][0] = std::complex<double>(0, -1.0);
            data[0][1] = -1.0;
            data[1][0] = std::complex<double>(0, -1.0);
            data[1][1] = 1.0;
            powOFsqrt2_inv = 1;
            break;
        case 112:
            data[0][0] = 1.0;
            data[0][1] = std::complex<double>(0, 1.0);
            data[1][0] = std::complex<double>(0, -1.0);
            data[1][1] = -1.0;
            powOFsqrt2_inv = 1;
            break;
        case 113:
            data[0][0] = std::complex<double>(-1.0, 1.0);
            data[0][1] = std::complex<double>(-1.0, 1.0);
            data[1][0] = std::complex<double>(1.0, 1.0);
            data[1][1] = std::complex<double>(-1.0, -1.0);
            powOFsqrt2_inv = 0;
            break;
        case 114:
            data[0][0] = 0.0;
            data[0][1] = std::complex<double>(1.0, 1.0);
            data[1][0] = std::complex<double>(-1.0, 1.0);
            data[1][1] = 0.0;
            powOFsqrt2_inv = 1;
            break;
        case 115:
            data[0][0] = std::complex<double>(1.0, -1.0);
            data[0][1] = std::complex<double>(1.0, -1.0);
            data[1][0] = std::complex<double>(-1.0, -1.0);
            data[1][1] = std::complex<double>(1.0, 1.0);
            powOFsqrt2_inv = 0;
            break;
        case 116:
            data[0][0] = std::complex<double>(-1.0, 1.0);
            data[0][1] = std::complex<double>(1.0, 1.0);
            data[1][0] = std::complex<double>(1.0, 1.0);
            data[1][1] = std::complex<double>(-1.0, 1.0);
            powOFsqrt2_inv = 0;
            break;
        case 117:
            data[0][0] = 0.0;
            data[0][1] = std::complex<double>(1.0, -1.0);
            data[1][0] = std::complex<double>(-1.0, -1.0);
            data[1][1] = 0.0;
            powOFsqrt2_inv = 1;
            break;
        case 118:
            data[0][0] = std::complex<double>(1.0, -1.0);
            data[0][1] = std::complex<double>(-1.0, -1.0);
            data[1][0] = std::complex<double>(-1.0, -1.0);
            data[1][1] = std::complex<double>(1.0, -1.0);
            powOFsqrt2_inv = 0;
            break;
        case 119:
            data[0][0] = std::complex<double>(1.0, 1.0);
            data[0][1] = 0.0;
            data[1][0] = 0.0;
            data[1][1] = std::complex<double>(-1.0, -1.0);
            powOFsqrt2_inv = 1;
            break;
        case 120:
            data[0][0] = 0.0;
            data[0][1] = std::complex<double>(-1.0, 1.0);
            data[1][0] = std::complex<double>(1.0, 1.0);
            data[1][1] = 0.0;
            powOFsqrt2_inv = 1;
            break;
        case 121:
            data[0][0] = std::complex<double>(-1.0, -1.0);
            data[0][1] = std::complex<double>(1.0, -1.0);
            data[1][0] = std::complex<double>(-1.0, 1.0);
            data[1][1] = std::complex<double>(1.0, 1.0);
            powOFsqrt2_inv = 0;
            break;
        case 122:
            data[0][0] = std::complex<double>(1.0, 1.0);
            data[0][1] = std::complex<double>(-1.0, -1.0);
            data[1][0] = std::complex<double>(-1.0, 1.0);
            data[1][1] = std::complex<double>(-1.0, 1.0);
            powOFsqrt2_inv = 0;
            break;
        case 123:
            data[0][0] = std::complex<double>(-1.0, -1.0);
            data[0][1] = 0.0;
            data[1][0] = 0.0;
            data[1][1] = std::complex<double>(-1.0, 1.0);
            powOFsqrt2_inv = 1;
            break;
        case 124:
            data[0][0] = std::complex<double>(-1.0, 1.0);
            data[0][1] = std::complex<double>(-1.0, -1.0);
            data[1][0] = std::complex<double>(-1.0, -1.0);
            data[1][1] = std::complex<double>(-1.0, 1.0);
            powOFsqrt2_inv = 0;
            break;
        case 125:
            data[0][0] = std::complex<double>(1.0, -1.0);
            data[0][1] = std::complex<double>(-1.0, 1.0);
            data[1][0] = std::complex<double>(-1.0, -1.0);
            data[1][1] = std::complex<double>(-1.0, -1.0);
            powOFsqrt2_inv = 0;
            break;
        case 126:
            data[0][0] = std::complex<double>(-1.0, 1.0);
            data[0][1] = std::complex<double>(-1.0, -1.0);
            data[1][0] = std::complex<double>(1.0, -1.0);
            data[1][1] = std::complex<double>(-1.0, -1.0);
            powOFsqrt2_inv = 0;
            break;
        case 127:
            data[0][0] = -1.0;
            data[0][1] = -1.0;
            data[1][0] = std::complex<double>(0, -1.0);
            data[1][1] = std::complex<double>(0, 1.0);
            powOFsqrt2_inv = 1;
            break;
        case 128:
            data[0][0] = 0.0;
            data[0][1] = std::complex<double>(-1.0, -1.0);
            data[1][0] = std::complex<double>(1.0, -1.0);
            data[1][1] = 0.0;
            powOFsqrt2_inv = 1;
            break;
        case 129:
            data[0][0] = std::complex<double>(-1.0, -1.0);
            data[0][1] = std::complex<double>(-1.0, 1.0);
            data[1][0] = std::complex<double>(1.0, -1.0);
            data[1][1] = std::complex<double>(1.0, 1.0);
            powOFsqrt2_inv = 0;
            break;
        case 130:
            data[0][0] = std::complex<double>(-1.0, 1.0);
            data[0][1] = 0.0;
            data[1][0] = 0.0;
            data[1][1] = std::complex<double>(-1.0, -1.0);
            powOFsqrt2_inv = 1;
            break;
        case 131:
            data[0][0] = 0.0;
            data[0][1] = std::complex<double>(-1.0, 1.0);
            data[1][0] = std::complex<double>(1.0, -1.0);
            data[1][1] = 0.0;
            powOFsqrt2_inv = 1;
            break;
        case 132:
            data[0][0] = std::complex<double>(1.0, -1.0);
            data[0][1] = 0.0;
            data[1][0] = 0.0;
            data[1][1] = std::complex<double>(1.0, 1.0);
            powOFsqrt2_inv = 1;
            break;
        case 133:
            data[0][0] = std::complex<double>(1.0, 1.0);
            data[0][1] = std::complex<double>(1.0, 1.0);
            data[1][0] = std::complex<double>(1.0, -1.0);
            data[1][1] = std::complex<double>(-1.0, 1.0);
            powOFsqrt2_inv = 0;
            break;
        case 134:
            data[0][0] = std::complex<double>(-1.0, 1.0);
            data[0][1] = std::complex<double>(-1.0, 1.0);
            data[1][0] = std::complex<double>(-1.0, -1.0);
            data[1][1] = std::complex<double>(1.0, 1.0);
            powOFsqrt2_inv = 0;
            break;
        case 135:
            data[0][0] = 0.0;
            data[0][1] = std::complex<double>(1.0, 1.0);
            data[1][0] = std::complex<double>(1.0, 1.0);
            data[1][1] = 0.0;
            powOFsqrt2_inv = 1;
            break;
        case 136:
            data[0][0] = std::complex<double>(1.0, -1.0);
            data[0][1] = std::complex<double>(1.0, -1.0);
            data[1][0] = std::complex<double>(1.0, 1.0);
            data[1][1] = std::complex<double>(-1.0, -1.0);
            powOFsqrt2_inv = 0;
            break;
        case 137:
            data[0][0] = 0.0;
            data[0][1] = std::complex<double>(1.0, -1.0);
            data[1][0] = std::complex<double>(-1.0, 1.0);
            data[1][1] = 0.0;
            powOFsqrt2_inv = 1;
            break;
        case 138:
            data[0][0] = std::complex<double>(1.0, 1.0);
            data[0][1] = 0.0;
            data[1][0] = 0.0;
            data[1][1] = std::complex<double>(-1.0, 1.0);
            powOFsqrt2_inv = 1;
            break;
        case 139:
            data[0][0] = 0.0;
            data[0][1] = std::complex<double>(-1.0, -1.0);
            data[1][0] = std::complex<double>(-1.0, -1.0);
            data[1][1] = 0.0;
            powOFsqrt2_inv = 1;
            break;
        case 140:
            data[0][0] = std::complex<double>(-1.0, -1.0);
            data[0][1] = std::complex<double>(1.0, 1.0);
            data[1][0] = std::complex<double>(1.0, -1.0);
            data[1][1] = std::complex<double>(1.0, -1.0);
            powOFsqrt2_inv = 0;
            break;
        case 141:
            data[0][0] = 0.0;
            data[0][1] = std::complex<double>(1.0, 1.0);
            data[1][0] = std::complex<double>(1.0, -1.0);
            data[1][1] = 0.0;
            powOFsqrt2_inv = 1;
            break;
        case 142:
            data[0][0] = std::complex<double>(-1.0, -1.0);
            data[0][1] = 0.0;
            data[1][0] = 0.0;
            data[1][1] = std::complex<double>(1.0, 1.0);
            powOFsqrt2_inv = 1;
            break;
        case 143:
            data[0][0] = std::complex<double>(-1.0, 1.0);
            data[0][1] = std::complex<double>(1.0, -1.0);
            data[1][0] = std::complex<double>(1.0, 1.0);
            data[1][1] = std::complex<double>(1.0, 1.0);
            powOFsqrt2_inv = 0;
            break;
        case 144:
            data[0][0] = std::complex<double>(0, 1.0);
            data[0][1] = 1.0;
            data[1][0] = 1.0;
            data[1][1] = std::complex<double>(0, 1.0);
            powOFsqrt2_inv = 1;
            break;
        case 145:
            data[0][0] = 0.0;
            data[0][1] = std::complex<double>(1.0, -1.0);
            data[1][0] = std::complex<double>(1.0, 1.0);
            data[1][1] = 0.0;
            powOFsqrt2_inv = 1;
            break;
        case 146:
            data[0][0] = std::complex<double>(-1.0, -1.0);
            data[0][1] = std::complex<double>(-1.0, -1.0);
            data[1][0] = std::complex<double>(-1.0, 1.0);
            data[1][1] = std::complex<double>(1.0, -1.0);
            powOFsqrt2_inv = 0;
            break;
        case 147:
            data[0][0] = std::complex<double>(-1.0, 1.0);
            data[0][1] = std::complex<double>(1.0, 1.0);
            data[1][0] = std::complex<double>(-1.0, 1.0);
            data[1][1] = std::complex<double>(-1.0, -1.0);
            powOFsqrt2_inv = 0;
            break;
        case 148:
            data[0][0] = std::complex<double>(-1.0, 1.0);
            data[0][1] = std::complex<double>(-1.0, 1.0);
            data[1][0] = std::complex<double>(1.0, -1.0);
            data[1][1] = std::complex<double>(-1.0, 1.0);
            powOFsqrt2_inv = 0;
            break;
        case 149:
            data[0][0] = std::complex<double>(1.0, -1.0);
            data[0][1] = std::complex<double>(-1.0, -1.0);
            data[1][0] = std::complex<double>(1.0, -1.0);
            data[1][1] = std::complex<double>(1.0, 1.0);
            powOFsqrt2_inv = 0;
            break;
        case 150:
            data[0][0] = std::complex<double>(0, 1.0);
            data[0][1] = -1.0;
            data[1][0] = -1.0;
            data[1][1] = std::complex<double>(0, 1.0);
            powOFsqrt2_inv = 1;
            break;
        case 151:
            data[0][0] = std::complex<double>(-1.0, -1.0);
            data[0][1] = std::complex<double>(1.0, 1.0);
            data[1][0] = std::complex<double>(1.0, 1.0);
            data[1][1] = std::complex<double>(1.0, 1.0);
            powOFsqrt2_inv = 0;
            break;
        case 152:
            data[0][0] = std::complex<double>(0, -1.0);
            data[0][1] = 1.0;
            data[1][0] = 1.0;
            data[1][1] = std::complex<double>(0, -1.0);
            powOFsqrt2_inv = 1;
            break;
        case 153:
            data[0][0] = std::complex<double>(1.0, -1.0);
            data[0][1] = std::complex<double>(1.0, -1.0);
            data[1][0] = std::complex<double>(-1.0, 1.0);
            data[1][1] = std::complex<double>(1.0, -1.0);
            powOFsqrt2_inv = 0;
            break;
        case 154:
            data[0][0] = std::complex<double>(1.0, 1.0);
            data[0][1] = std::complex<double>(1.0, -1.0);
            data[1][0] = std::complex<double>(1.0, 1.0);
            data[1][1] = std::complex<double>(-1.0, 1.0);
            powOFsqrt2_inv = 0;
            break;
        case 155:
            data[0][0] = std::complex<double>(1.0, 1.0);
            data[0][1] = std::complex<double>(-1.0, -1.0);
            data[1][0] = std::complex<double>(-1.0, -1.0);
            data[1][1] = std::complex<double>(-1.0, -1.0);
            powOFsqrt2_inv = 0;
            break;
        case 156:
            data[0][0] = -1.0;
            data[0][1] = std::complex<double>(0, 1.0);
            data[1][0] = std::complex<double>(0, -1.0);
            data[1][1] = 1.0;
            powOFsqrt2_inv = 1;
            break;
        case 157:
            data[0][0] = std::complex<double>(-1.0, 1.0);
            data[0][1] = std::complex<double>(1.0, 1.0);
            data[1][0] = std::complex<double>(1.0, -1.0);
            data[1][1] = std::complex<double>(1.0, 1.0);
            powOFsqrt2_inv = 0;
            break;
        case 158:
            data[0][0] = std::complex<double>(-1.0, -1.0);
            data[0][1] = std::complex<double>(-1.0, -1.0);
            data[1][0] = std::complex<double>(-1.0, -1.0);
            data[1][1] = std::complex<double>(1.0, 1.0);
            powOFsqrt2_inv = 0;
            break;
        case 159:
            data[0][0] = -1.0;
            data[0][1] = std::complex<double>(0, -1.0);
            data[1][0] = std::complex<double>(0, 1.0);
            data[1][1] = 1.0;
            powOFsqrt2_inv = 1;
            break;
        case 160:
            data[0][0] = std::complex<double>(-1.0, -1.0);
            data[0][1] = std::complex<double>(1.0, -1.0);
            data[1][0] = std::complex<double>(1.0, 1.0);
            data[1][1] = std::complex<double>(1.0, -1.0);
            powOFsqrt2_inv = 0;
            break;
        case 161:
            data[0][0] = std::complex<double>(0, -1.0);
            data[0][1] = -1.0;
            data[1][0] = -1.0;
            data[1][1] = std::complex<double>(0, -1.0);
            powOFsqrt2_inv = 1;
            break;
        case 162:
            data[0][0] = std::complex<double>(-1.0, 1.0);
            data[0][1] = std::complex<double>(1.0, -1.0);
            data[1][0] = std::complex<double>(-1.0, 1.0);
            data[1][1] = std::complex<double>(-1.0, 1.0);
            powOFsqrt2_inv = 0;
            break;
        case 163:
            data[0][0] = std::complex<double>(1.0, -1.0);
            data[0][1] = std::complex<double>(-1.0, -1.0);
            data[1][0] = std::complex<double>(-1.0, 1.0);
            data[1][1] = std::complex<double>(-1.0, -1.0);
            powOFsqrt2_inv = 0;
            break;
        case 164:
            data[0][0] = std::complex<double>(1.0, -1.0);
            data[0][1] = std::complex<double>(-1.0, 1.0);
            data[1][0] = std::complex<double>(1.0, -1.0);
            data[1][1] = std::complex<double>(1.0, -1.0);
            powOFsqrt2_inv = 0;
            break;
        case 165:
            data[0][0] = std::complex<double>(-1.0, 1.0);
            data[0][1] = std::complex<double>(-1.0, -1.0);
            data[1][0] = std::complex<double>(-1.0, 1.0);
            data[1][1] = std::complex<double>(1.0, 1.0);
            powOFsqrt2_inv = 0;
            break;
        case 166:
            data[0][0] = std::complex<double>(1.0, 1.0);
            data[0][1] = std::complex<double>(-1.0, 1.0);
            data[1][0] = std::complex<double>(-1.0, -1.0);
            data[1][1] = std::complex<double>(-1.0, 1.0);
            powOFsqrt2_inv = 0;
            break;
        case 167:
            data[0][0] = std::complex<double>(1.0, -1.0);
            data[0][1] = std::complex<double>(1.0, 1.0);
            data[1][0] = std::complex<double>(1.0, -1.0);
            data[1][1] = std::complex<double>(-1.0, -1.0);
            powOFsqrt2_inv = 0;
            break;
        case 168:
            data[0][0] = std::complex<double>(1.0, 1.0);
            data[0][1] = std::complex<double>(-1.0, -1.0);
            data[1][0] = std::complex<double>(1.0, 1.0);
            data[1][1] = std::complex<double>(1.0, 1.0);
            powOFsqrt2_inv = 0;
            break;
        case 169:
            data[0][0] = std::complex<double>(1.0, -1.0);
            data[0][1] = std::complex<double>(-1.0, 1.0);
            data[1][0] = std::complex<double>(-1.0, 1.0);
            data[1][1] = std::complex<double>(-1.0, 1.0);
            powOFsqrt2_inv = 0;
            break;
        case 170:
            data[0][0] = std::complex<double>(-1.0, -1.0);
            data[0][1] = std::complex<double>(-1.0, 1.0);
            data[1][0] = std::complex<double>(-1.0, -1.0);
            data[1][1] = std::complex<double>(1.0, -1.0);
            powOFsqrt2_inv = 0;
            break;
        case 171:
            data[0][0] = std::complex<double>(1.0, 1.0);
            data[0][1] = std::complex<double>(1.0, 1.0);
            data[1][0] = std::complex<double>(-1.0, -1.0);
            data[1][1] = std::complex<double>(1.0, 1.0);
            powOFsqrt2_inv = 0;
            break;
        case 172:
            data[0][0] = std::complex<double>(-1.0, 1.0);
            data[0][1] = 0.0;
            data[1][0] = 0.0;
            data[1][1] = std::complex<double>(-1.0, 1.0);
            powOFsqrt2_inv = 1;
            break;
        case 173:
            data[0][0] = 0.0;
            data[0][1] = std::complex<double>(-1.0, -1.0);
            data[1][0] = std::complex<double>(-1.0, 1.0);
            data[1][1] = 0.0;
            powOFsqrt2_inv = 1;
            break;
        case 174:
            data[0][0] = std::complex<double>(1.0, -1.0);
            data[0][1] = 0.0;
            data[1][0] = 0.0;
            data[1][1] = std::complex<double>(1.0, -1.0);
            powOFsqrt2_inv = 1;
            break;
        case 175:
            data[0][0] = std::complex<double>(-1.0, 1.0);
            data[0][1] = 0.0;
            data[1][0] = 0.0;
            data[1][1] = std::complex<double>(1.0, 1.0);
            powOFsqrt2_inv = 1;
            break;
        case 176:
            data[0][0] = 0.0;
            data[0][1] = std::complex<double>(-1.0, 1.0);
            data[1][0] = std::complex<double>(-1.0, -1.0);
            data[1][1] = 0.0;
            powOFsqrt2_inv = 1;
            break;
        case 177:
            data[0][0] = std::complex<double>(1.0, -1.0);
            data[0][1] = 0.0;
            data[1][0] = 0.0;
            data[1][1] = std::complex<double>(-1.0, -1.0);
            powOFsqrt2_inv = 1;
            break;
        case 178:
            data[0][0] = std::complex<double>(1.0, 1.0);
            data[0][1] = 0.0;
            data[1][0] = 0.0;
            data[1][1] = std::complex<double>(1.0, 1.0);
            powOFsqrt2_inv = 1;
            break;
        case 179:
            data[0][0] = 0.0;
            data[0][1] = std::complex<double>(-1.0, 1.0);
            data[1][0] = std::complex<double>(-1.0, 1.0);
            data[1][1] = 0.0;
            powOFsqrt2_inv = 1;
            break;
        case 180:
            data[0][0] = std::complex<double>(-1.0, -1.0);
            data[0][1] = 0.0;
            data[1][0] = 0.0;
            data[1][1] = std::complex<double>(1.0, -1.0);
            powOFsqrt2_inv = 1;
            break;
        case 181:
            data[0][0] = 0.0;
            data[0][1] = std::complex<double>(1.0, 1.0);
            data[1][0] = std::complex<double>(-1.0, -1.0);
            data[1][1] = 0.0;
            powOFsqrt2_inv = 1;
            break;
        case 182:
            data[0][0] = 0.0;
            data[0][1] = std::complex<double>(1.0, -1.0);
            data[1][0] = std::complex<double>(1.0, -1.0);
            data[1][1] = 0.0;
            powOFsqrt2_inv = 1;
            break;
        case 183:
            data[0][0] = std::complex<double>(-1.0, 1.0);
            data[0][1] = 0.0;
            data[1][0] = 0.0;
            data[1][1] = std::complex<double>(1.0, -1.0);
            powOFsqrt2_inv = 1;
            break;
        case 184:
            data[0][0] = 0.0;
            data[0][1] = std::complex<double>(-1.0, -1.0);
            data[1][0] = std::complex<double>(1.0, 1.0);
            data[1][1] = 0.0;
            powOFsqrt2_inv = 1;
            break;
        case 185:
            data[0][0] = std::complex<double>(1.0, -1.0);
            data[0][1] = 0.0;
            data[1][0] = 0.0;
            data[1][1] = std::complex<double>(-1.0, 1.0);
            powOFsqrt2_inv = 1;
            break;
        case 186:
            data[0][0] = std::complex<double>(-1.0, 1.0);
            data[0][1] = std::complex<double>(1.0, -1.0);
            data[1][0] = std::complex<double>(1.0, -1.0);
            data[1][1] = std::complex<double>(1.0, -1.0);
            powOFsqrt2_inv = 0;
            break;
        case 187:
            data[0][0] = std::complex<double>(-1.0, -1.0);
            data[0][1] = 0.0;
            data[1][0] = 0.0;
            data[1][1] = std::complex<double>(-1.0, -1.0);
            powOFsqrt2_inv = 1;
            break;
        case 188:
            data[0][0] = std::complex<double>(-1.0, 1.0);
            data[0][1] = std::complex<double>(-1.0, 1.0);
            data[1][0] = std::complex<double>(-1.0, 1.0);
            data[1][1] = std::complex<double>(1.0, -1.0);
            powOFsqrt2_inv = 0;
            break;
        case 189:
            data[0][0] = std::complex<double>(-1.0, -1.0);
            data[0][1] = std::complex<double>(-1.0, -1.0);
            data[1][0] = std::complex<double>(1.0, 1.0);
            data[1][1] = std::complex<double>(-1.0, -1.0);
            powOFsqrt2_inv = 0;
            break;
        case 190:
            data[0][0] = std::complex<double>(1.0, -1.0);
            data[0][1] = std::complex<double>(1.0, -1.0);
            data[1][0] = std::complex<double>(1.0, -1.0);
            data[1][1] = std::complex<double>(-1.0, 1.0);
            powOFsqrt2_inv = 0;
            break;
        case 191:
            data[0][0] = std::complex<double>(-1.0, -1.0);
            data[0][1] = std::complex<double>(1.0, 1.0);
            data[1][0] = std::complex<double>(-1.0, -1.0);
            data[1][1] = std::complex<double>(-1.0, -1.0);
            powOFsqrt2_inv = 0;
            break;
        default:
            throw std::runtime_error("Matrix index out of range");
        }
    }
