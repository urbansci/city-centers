# -*- coding:utf-8 -*-
# @author  : Shuai Pang
# @time    : 2025-11

"""
Module storages model parameters and predefined terms, etc.
"""

# File path
FILE_PATH = "D:\\Research\\CityCenter\\"

# Method parameter
MINIMUM_CONTOUR_AREA = None  # in units of km

# Income group by World Bank for 2020
GROUPS = [ 'Low income', 'Lower-middle income', 'Upper-middle income','High income',]
GROUPS_DICT = {'High income': [0, 1, 2, 4, 5, 6, 9, 11, 12, 14, 17, 18, 29, 31, 34, 38, 40, 45, 50, 55, 59, 65, 107, 133, 134, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 150, 151, 152, 153, 154, 155, 156, 157, 158, 160, 161, 162, 163, 164, 169, 170, 171, 172, 174, 179, 180, 181, 182, 183, 191, 196, 201, 202, 204, 205, 211, 217, 219, 220, 221, 226],
               'Low income': [72, 74, 75, 78, 79, 80, 81, 82, 84, 85, 88, 90, 92, 93, 95, 100, 101, 103, 104, 108, 111, 112, 187, 189, 209, 230],
               'Lower-middle income': [10, 20, 25, 27, 28, 39, 41, 42, 46, 61, 62, 66, 67, 69, 70, 73, 76, 77, 83, 86, 87, 89, 91, 94, 97, 98, 102, 105, 106, 110, 114, 118, 119, 122, 129, 167, 185, 192, 193, 194, 197, 199, 206, 207, 210, 213, 214, 216, 218, 223, 224, 225, 228, 231, 233],
               'Upper-middle income': [7, 8, 13, 16, 22, 23, 26, 30, 32, 35, 36, 37, 44, 47, 48, 57, 68, 96, 99, 113, 116, 117, 120, 124, 125, 127, 128, 130, 131, 132, 136, 138, 149, 159, 165, 166, 168, 173, 175, 176, 178, 184, 186, 188, 190, 195, 200, 203, 208, 212, 222, 227, 229, 232]}
GROUPS_MAP = {idx: group for group in GROUPS_DICT.keys() for idx in GROUPS_DICT[group]}

# Continent partition
CONTS = ['SA', 'OC', 'NA', 'EU', 'AS', 'AF']
CONTS_DICT = {'NA': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36],
 'OC': [37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65],
 'AF': [66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122],
 'SA': [123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136],
 'EU': [137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183],
 'AS': [184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233]
}
CONTS_MAP = {idx: cont for cont in CONTS_DICT.keys() for idx in CONTS_DICT[cont]}

# Region partition
REGIONS = ['Sub-Saharan Africa', 'South Asia', 'North America', 'M. East & N. Africa', 'Latin Am. & Carib.', 'Europe & Central Asia', 'East Asia & Pacific']
REGIONS_DICT = {'East Asia & Pacific': [44, 45, 196, 59, 219, 220, 40, 50, 65, 55, 38, 211, 230, 213, 39, 194, 62, 218, 42, 225, 61, 41, 197, 46, 231, 232, 47, 224, 222, 37, 199, 227, 48, 57],
                'Europe & Central Asia': [174, 160, 154, 173, 169, 205, 157, 148, 145, 142, 143, 155, 150, 180, 179, 0, 161, 140, 153, 172, 146, 163, 147, 156, 171, 151, 141, 152, 182, 165, 138, 170, 158, 164, 183, 139, 162, 144, 185, 216, 228, 178, 186, 208, 149, 168, 203, 212, 159, 176, 175, 166, 229, 184, 167],
                'Latin Am. & Carib.': [18, 31, 5, 29, 12, 9, 134, 124, 36, 11, 17, 34, 6, 133, 14, 129, 10, 28, 25, 132, 20, 128, 136, 35, 7, 23, 13, 127, 27, 32, 22, 16, 8, 131, 130, 26, 30, 125, 123],
                'M. East & N. Africa': [201, 204, 191, 181, 226, 202, 217, 221, 187, 209, 83, 70, 200, 188, 69, 67, 233, 66, 206, 190, 68],
                'North America': [4, 1, 2],
                'South Asia': [189, 210, 193, 223, 192, 207, 214, 195],
                'Sub-Saharan Africa': [107, 80, 104, 88, 85, 103, 81, 71, 78, 79, 82, 93, 111, 108, 72, 112, 74, 101, 90, 100, 75, 92, 95, 106, 91, 77, 94, 122, 97, 87, 118, 89, 84, 102, 119, 73, 86, 98, 76, 105, 110, 114, 117, 96, 99, 113, 116, 120]}
REGIONS_MAP = {idx: region for region in REGIONS_DICT.keys() for idx in REGIONS_DICT[region]}

# Socioeconomic indicator
COLUMNS = ['Count', 'Area', 'Pop', 'Bright']
COLUMN_NAMES = ['Count', 'Area', 'Population', 'Brightness']

# City categories
LABELS = ['mono', 'low poly', 'moderate poly', 'high poly']

# Country names
NAMES = {-1: 'World', 0: 'Greenland', 1: 'Canada', 2: 'US', 3: 'Saint Pierre And Miquelon', 4: 'Bermuda', 5: 'Bahamas', 6: 'Turks And Caicos Islands', 7: 'Cuba', 8: 'Mexico', 9: 'Cayman Islands', 10: 'Haiti', 11: 'Puerto Rico', 12: 'Virgin Islands,British', 13: 'Dominican Republic', 14: 'Virgin Islands,U.S.', 15: 'Anguilla', 16: 'Jamaica', 17: 'Saint Kitts And Nevis', 18: 'Antigua And Barbuda', 19: 'Montserrat', 20: 'Belize', 21: 'Guadeloupe', 22: 'Guatemala', 23: 'Dominica', 24: 'Martinique', 25: 'Nicaragua', 26: 'Saint Lucia', 27: 'El Salvador', 28: 'Honduras', 29: 'Barbados', 30: 'Saint Vincent And The Grenadines', 31: 'Aruba', 32: 'Grenada', 33: 'Netherlands Antilles', 34: 'Trinidad And Tobago', 35: 'Costa Rica', 36: 'Panama', 37: 'Marshall Islands', 38: 'Palau', 39: 'Kiribati', 40: 'Nauru', 41: 'Solomon Islands', 42: 'Papua New Guinea', 43: 'Cook Islands', 44: 'American Samoa', 45: 'Australia', 46: 'Vanuatu', 47: 'Fiji', 48: 'Tonga', 49: 'Niue', 50: 'New Caledonia', 51: 'Pitcairn', 52: 'Norfolk Island', 53: 'Heard Island And Mcdonald Islands', 54: 'Bouvet Island', 55: 'Northern Mariana Islands', 56: 'South Georgia And The South Sandwich Islands', 57: 'Tuvalu', 58: 'Tokelau', 59: 'French Polynesia', 60: 'French Southern Territories', 61: 'Samoa', 62: 'Federated States of Micronesia', 63: 'Wallis And Futuna', 64: 'Cocos(Keeling) Islands', 65: 'New Zealand', 66: 'Algeria', 67: 'Tunisia', 68: 'Libya', 69: 'Morocco', 70: 'Egypt', 71: 'Western Sahara', 72: 'Mali', 73: 'Mauritania', 74: 'Niger', 75: 'Sudan', 76: 'Senegal', 77: 'Cape Verde', 78: 'Ethiopia', 79: 'Gambia', 80: 'Burkina Faso', 81: 'Eritrea', 82: 'Guinea-Bissau', 83: 'Djibouti', 84: 'Guinea', 85: 'Chad', 86: 'Nigeria', 87: 'Ivory Coast', 88: 'Central African Republic', 89: 'Ghana', 90: 'Sierra Leone', 91: 'Benin', 92: 'Togo', 93: 'Liberia', 94: 'Cameroon', 95: 'Uganda', 96: 'Equatorial Guinea', 97: 'Congo', 98: 'Sao Tome And Principe', 99: 'Gabon', 100: 'Somalia', 101: 'Rwanda', 102: 'Kenya', 103: 'Democratic Republic of the Congo', 104: 'Burundi', 105: 'Tanzania', 106: 'Angola', 107: 'Seychelles', 108: 'Malawi', 109: 'Mayotte', 110: 'Zambia', 111: 'Madagascar', 112: 'Mozambique', 113: 'Mauritius', 114: 'Zimbabwe', 115: 'Reunion', 116: 'Namibia', 117: 'Botswana', 118: 'Swaziland', 119: 'Lesotho', 120: 'South Africa', 121: 'Saint Helena', 122: 'Comoros', 123: 'Venezuela', 124: 'Guyana', 125: 'Suriname', 126: 'French Guiana', 127: 'Ecuador', 128: 'Brazil', 129: 'Bolivia', 130: 'Peru', 131: 'Paraguay', 132: 'Argentina', 133: 'Uruguay', 134: 'Chile', 135: 'Falkland Islands(Malvinas)', 136: 'Colombia', 137: 'Svalbard And Jan Mayen', 138: 'Russia', 139: 'Sweden', 140: 'Iceland', 141: 'Norway', 142: 'Faroe Islands', 143: 'Finland', 144: 'UK', 145: 'Estonia', 146: 'Latvia', 147: 'Lithuania', 148: 'Denmark', 149: 'Belarus', 150: 'Germany', 151: 'Netherlands', 152: 'Poland', 153: 'Ireland', 154: 'Belgium', 155: 'France', 156: 'Luxembourg', 157: 'Czech Republic', 158: 'Slovakia', 159: 'Moldova', 160: 'Austria', 161: 'Hungary', 162: 'Switzerland', 163: 'Liechtenstein', 164: 'Slovenia', 165: 'Romania', 166: 'Serbia', 167: 'Ukraine', 168: 'Bosnia And Herzegovina', 169: 'Croatia', 170: 'San Marino', 171: 'Monaco', 172: 'Italy', 173: 'Bulgaria', 174: 'Andorra', 175: 'North Macedonia', 176: 'Montenegro', 177: 'Vatican', 178: 'Albania', 179: 'Greece', 180: 'Gibraltar', 181: 'Malta', 182: 'Portugal', 183: 'Spain', 184: 'Turkmenistan', 185: 'Kyrgyzstan', 186: 'Armenia', 187: 'Syrian Arab Republic', 188: 'Lebanon', 189: 'Afghanistan', 190: 'Iraq', 191: 'Kuwait', 192: 'Nepal', 193: 'Bhutan', 194: 'Laos', 195: 'Maldives', 196: 'Brunei', 197: 'East Timor', 198: 'Christmas Island', 199: 'Mongolia', 200: 'Jordan', 201: 'Bahrain', 202: 'Qatar', 203: 'Georgia', 204: 'Israel', 205: 'Cyprus', 206: 'Iran', 207: 'Pakistan', 208: 'Azerbaijan', 209: 'Yemen', 210: 'Bangladesh', 211: 'Singapore', 212: 'Kazakhstan', 213: 'Cambodia', 214: 'Srilanka', 215: 'British Indian Ocean Territory', 216: 'Tajikistan', 217: 'Saudi Arabia', 218: 'Myanmar', 219: 'Japan', 220: 'South Korea', 221: 'United Arab Emirates', 222: 'Malaysia', 223: 'India', 224: 'Indonesia', 225: 'Philippines', 226: 'Oman', 227: 'Thailand', 228: 'Uzbekistan', 229: 'Turkey', 230: "North Korea", 231: 'Viet Nam', 232: 'China', 233: 'Palestine'}

# RTree
contour_tree = None
poi_tree = None