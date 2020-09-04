import argparse
import json
import requests
from requests.auth import HTTPBasicAuth


def _get_connection(args):
    authobj = None
    if args.username is not None:
        authobj = HTTPBasicAuth(args.username, args.password)
    url = "%s/tradition/%s/section/%s" % (args.repository, args.tradition, args.section)
    return url, authobj


def get_baselist(args):
    if args.local:
        # Read the list in from the file specified
        with open(args.base, encoding="utf-8") as f:
            baselist = json.load(f)
    else:
        # Fetch the list from the server
        (url, authobj) = _get_connection(args)
        payload = {'final': 'true'}
        r = requests.get("%s/lemmareadings", params=payload, auth=authobj)
        r.raise_for_status()
        baselist = r.json()
    return baselist


def get_variants(args):
    if args.local:
        # Read the list in from the file specified
        with open(args.variantlist, encoding="utf-8") as f:
            variantlist = json.load(f)
    else:
        # Get the server URL
        (url, authobj) = _get_connection(args)
        # Get the relevant query parameters to pass
        payload = {}
        for ostr in args.option:
            opt = ostr.split('=')
            payload[opt[0]] = opt[1]
        r = requests.get("%s/lemmareadings", params=payload, auth=authobj)
        r.raise_for_status()
        variantlist = r.json()
    return variantlist


def make_applookup(variantlist):
    """Make a lookup table of lemma reading ID -> list of variants that need to be footnoted at that lemma."""
    applookup = {}
    for v in variantlist.get('variantlist'):
        # If the base chain exists, hash on the last reading in the base chain
        if 'base' in v:
            rid = v.get('base')[-1].get('id')
        # Otherwise hash on the end reading
        else:
            rid = v.get('after').get('id')
        # Add this variant to the list, or create a new list with the variant
        if rid in applookup:
            applookup.get(rid).append(v)
        else:
            applookup[rid] = [v]
    return applookup


def _is_ascii(s):
    for c in s:
        if ord(c) > 128:
            return False
    return True


def _lwrap(langtag, *rdgs):
    st = ''
    join_next = False
    for r in rdgs:
        rtext = r.get('normal_form', r.get('text'))
        if join_next or r.get('join_prior', False) or len(st) == 0:
            st += rtext
        else:
            st += ' ' + rtext
        join_next = r.get('join_next', False)
    # Cheap hack - is the string entirely ASCII?

    return st if (langtag is None or _is_ascii(st)) else "\\%s{%s}" % (langtag, st)


def _get_lemmareadings(vloc):
    """Returns an appropriate lemma text based on the variant location features."""
    if 'base' in vloc:
        # If there is a list of base readings, use it.
        # (Here is where we can shorten lemma strings if we see fit
        return vloc.get('base')
    else:
        # Otherwise this is an addition, so use the before and after readings.
        return [vloc.get('before'), vloc.get('after')]


def _make_witstr(witlist):
    sigla = []
    for k in witlist:
        if k == 'witnesses':
            sigla.extend(witlist[k])
        else:
            sigla.extend(["%s (%s)" % (x, k) for x in witlist[k]])
    return ' '.join(sorted(sigla))


def _get_vreadings(vloc, langtag, special):
    """Returns a string listing all the readings and their witnesses in a given variant location."""
    stringified = []
    for v in vloc.get('variants'):
        wits = _make_witstr(v.get('witnesses'))
        if 'readings' not in v:
            txt = '\\emph{om.}'
        elif special == 'transp' and v.get('displaced', False):
            # Is it pre or post?
            txt = '\\emph{transp. '
            anchorrank = v.get('anchor').get('rank')
            ourrank = v.get('readings')[0].get('rank')
            txt += 'prae} ' if ourrank < anchorrank else 'post} '
            txt += _lwrap(langtag, v.get('anchor'))
        else:
            txt = _lwrap(langtag, *v.get('readings'))
        if special == 'add':
            txt += ' \\emph{add.}'
        stringified.append("%s %s" % (wits, txt))
    return '; '.join(stringified)


def generate_latex(baselist, applist, langtag):
    output = ''
    join_next = False
    for r in baselist:
        if not (join_next or r.get('join_prior', False) or len(output) == 0):
            output += ' '
        if r.get('id') in applist:
            vhere = applist.get(r.get('id'))
            # Generate the apparatus entry.
            output += "\\edtext{%s}{" % _lwrap(langtag, r)
            for vloc in vhere:
                if 'base' not in vloc:
                    special = 'add'
                elif vloc.get('has_displacement', False):
                    special = 'transp'
                else:
                    special = None
                output += "{\\lemma{%s} \\Afootnote{%s}}" % (
                    _lwrap(langtag, *_get_lemmareadings(vloc)), _get_vreadings(vloc, langtag, special))
            output += '}'
        else:
            output += _lwrap(langtag, r)
        join_next = r.get('join_next', False)
    print(output)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--remote",
        action="store_true",
        help="Get the data from a remote URL"
    )
    parser.add_argument(
        "--local",
        action="store_true",
        help="Get the data from local JSON files"
    )

    remote = parser.add_argument_group('Remote server connection')
    remote.add_argument(
        "-r",
        "--repository",
        help="URL to tradition repository"
    )
    remote.add_argument(
        "-u",
        "--username",
        help="HTTP basic auth username for tradition repository"
    )
    remote.add_argument(
        "-p",
        "--password",
        help="HTTP basic auth password for tradition repository"
    )
    remote.add_argument(
        "-t",
        "--tradition-id",
        help="ID of tradition to process"
    )
    remote.add_argument(
        "-s",
        "--section",
        help="Restrict processing to given section"
    )
    remote.add_argument(
        "-o",
        "--option",
        action="append",
        help="Set a query parameter ('name=value') to generate the variant list"
    )

    local = parser.add_argument_group('Data from local files')
    local.add_argument(
        "-b",
        "--base",
        help="Path to file containing a list of base readings in Stemmarest JSON format"
    )
    local.add_argument(
        "-v",
        "--variantlist",
        help="Path to file containing a list of variants in Stemmarest JSON format"
    )

    parser.add_argument(
        "--language",
        help="Polyglossia language tag to use for wrapping readings"
    )

    args = parser.parse_args()
    # Sanity check the arguments - we either need both files for local run, or a full URL spec for remote run
    if args.remote and args.local:
        raise Exception("Cannot specify both local and remote runs")
    if args.remote:
        # We need at least a repo URL, a tradition ID, and a section ID.
        if args.repository is None or args.tradition is None or args.section is None:
            raise Exception("For remote runs please specify repository URL and tradition / section IDs")
    else:
        if args.base is None or args.variantlist is None:
            raise Exception("For local runs please specify a file containing base readings, "
                            "and a file containing variants")

    # Get the list of base readings
    baselist = get_baselist(args)
    # Read in the variants and hash them accordingly
    applist = make_applookup(get_variants(args))
    # Spit out the LaTeX lines
    generate_latex(baselist, applist, args.language)
