import pymel.core as pm

class channelController:
    def create_gradient_params(self, channel, node_name):
        self.addControl('%sChannelMode' % channel, label='Channel Mode')

        self.beginLayout('Channel Controls', collapse=True)

        self.create_float_ramp(node_name, channel)
        self.addControl('%sContrast' % channel, label='Contrast')
        self.addControl('%sContrastPivot' % channel, label='Contrast Pivot')
        self.addControl('%sBias' % channel, label='Bias')
        self.addControl('%sGain' % channel, label='Gain')
        self.addControl('%sOutputMin' % channel, label='Output Min')
        self.addControl('%sOutputMax' % channel, label='Output Max')
        self.addControl('%sClampMin' % channel, label='Clamp Min')
        self.addControl('%sClampMax' % channel, label='Clamp Max')

        self.addSeparator()
        self.create_rgb_ramp(node_name, channel)
        self.addControl('%sGamma' % channel, label='Gamma')
        self.addControl('%sHueShift' % channel, label='Hue Shift')
        self.addControl('%sSaturation' % channel, label='Saturation')
        self.addControl('%sContrast' % channel, label='Contrast')
        self.addControl('%sContrastPivot' % channel, label='Contrast Pivot')
        self.addControl('%sExposure' % channel, label='Exposure')
        self.addControl('%sMultiply' % channel, label='Multiply')
        self.addControl('%sAdd' % channel, label='Add')

        self.endLayout()

    def create_float_ramp(self, channel, node_name):
        pm.mel.eval('source AEaddRampControl.mel; AEaddRampControl("%s.%sFloatRamp")' % (node_name, channel))
        self.addControl('%sInputMin' % channel, label='Ramp Input Min')
        self.addControl('%sInputMax' % channel, label='Ramp Input Max')

    def create_rgb_ramp(self, channel, node_name):
        pm.mel.eval('source AEaddRampControl.mel; AEaddRampControl("%s.%sRgbRamp")' % (node_name, channel))
        self.addControl('%sInputMin' % channel, label='Ramp Input Min')
        self.addControl('%sInputMax' % channel, label='Ramp Input Max')
