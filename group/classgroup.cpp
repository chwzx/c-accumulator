//
// Created by root on 23-4-18.
//

#include "classgroup.h"
#include "../util.h"

classgroup::classElem classgroup::id() {
    return classgroup::id_(DETAL);
}

classgroup::classElem classgroup::op(classgroup::classElem a, classgroup::classElem b) {
    return classgroup::op_(a, b);
}

classgroup::classElem classgroup::exp(classgroup::classElem a, mpz_class n) {
    return classgroup::exp_(a, n);
}

classgroup::classElem classgroup::inv(classgroup::classElem a) {
    return classgroup::inv_(a);
}

classgroup::classElem classgroup::unknown_order_elem() {
    return classgroup::unknown_order_elem_(DETAL);
}

classgroup::classElem classgroup::normalize(mpz_class a, mpz_class b, mpz_class c) {
    if(is_normal(a, b, c)){
        return {a, b, c};
    }
    // r = floor_div((a - b), 2a)
    // (a, b, c) = (a, b + 2ra, ar^2 + br + c)
    mpz_class r, rem, s = a - b, t = 2 * a;
    mpz_div(r.get_mpz_t(), s.get_mpz_t(), t.get_mpz_t());
    mpz_class new_b = b + 2 * (r * a);
    mpz_class new_c = c + b * r + a * r * r;
    return {a, new_b, new_c};
}

classgroup::classElem classgroup::reduce(mpz_class a, mpz_class b, mpz_class c) {
    while(!is_reduced(a, b, c)){
        // s = floor_div(c + b, 2c)
        mpz_class m = c + b, n = 2 * c, s;
        mpz_div(s.get_mpz_t(), m.get_mpz_t(), n.get_mpz_t());
        // (a, b, c) = (c, −b + 2sc, cs^2 − bs + a)
        mpz_class old_a = a;
        mpz_class old_b = b;
        a = c;
        b = -b + 2 * (s * c);
        c = -(old_b * s) + old_a + c * s * s;
    }
    return normalize(a, b, c);
}

classgroup::classElem classgroup::square(classgroup::classElem x) {
    util u;
    mpz_class mu = u.solve_linear_congruence(x.b, x.c, x.a).value().first;

    // A = a^2
    // B = b - 2a * mu
    // tmp = (b * mu - c) / a
    // C = mu^2 - tmp
    mpz_class a = x.a * x.a;
    mpz_class b = x.b - 2 * x.a * mu;
    mpz_class tmp = (x.b * mu - x.c) / x.a;
    mpz_class c = mu * mu - tmp;

    return elem({a, b, c});
}

mpz_class classgroup::discriminant(mpz_class a, mpz_class b, mpz_class c) {
    return b * b - 4 * a * c;
}

bool classgroup::validate(mpz_class a, mpz_class b, mpz_class c) {
    return discriminant(a, b, c) == DETAL;
}

bool classgroup::is_reduced(mpz_class a, mpz_class b, mpz_class c) {
    return is_normal(a, b, c) && (a <=c && !(a == c && b < 0));
}

bool classgroup::is_normal(mpz_class a, mpz_class b, mpz_class c) {
    mpz_class neg_a;
    neg_a = -a;
    if(neg_a < b  && b <= a)
        return true;
    return false;
}

classgroup::classElem classgroup::op_(classgroup::classElem x, classgroup::classElem y) {
    // g = (b1 + b2) / 2
    // h = (b2 - b1) / 2
    // w = gcd(a1, a2, g)
    mpz_class g = x.b + y.b;
    mpz_div_ui(g.get_mpz_t(), g.get_mpz_t(), 2);
    mpz_class h = y.b - x.b;
    mpz_div_ui(h.get_mpz_t(), h.get_mpz_t(), 2);
    mpz_class w = 1;
    mpz_gcd(w.get_mpz_t(), x.a.get_mpz_t(), y.a.get_mpz_t());
    mpz_gcd(w.get_mpz_t(), w.get_mpz_t(), g.get_mpz_t());

    // j = w
    // s = a1 / w
    // t = a2 / w
    // u = g / w
    // r = 0
    mpz_class j = w;
    mpz_class s = 1;
    mpz_div(s.get_mpz_t(), x.a.get_mpz_t(), w.get_mpz_t());
    mpz_class t = 1;
    mpz_div(t.get_mpz_t(), y.a.get_mpz_t(), w.get_mpz_t());
    mpz_class u = 1;
    mpz_div(u.get_mpz_t(), g.get_mpz_t(), w.get_mpz_t());

    // a = tu
    // b = hu + sc
    // m = st
    // Solve linear congruence `(tu)k = hu + sc mod st` or `ak = b mod m` for solutions `k`.
    mpz_class a = t * u;
    mpz_class b = h * u + (s * x.c);
    mpz_class m = s * t;
    util ul;
    mpz_class mu = ul.solve_linear_congruence(a, b, m).value().first;
    mpz_class v = ul.solve_linear_congruence(a, b, m).value().second;
    // a = tv
    // b = h - t * mu
    // m = s
    // Solve linear congruence `(tv)k = h - t * mu mod s` or `ak = b mod m` for solutions `k`.
    a = t * v;
    b = h - t * mu;
    m = s;
    mpz_class lambda = ul.solve_linear_congruence(a, b, m).value().first;;
    // k = mu + v * lambda
    // l = (k * t - h) / s
    // m = (tuk - hu - cs) / st
    mpz_class k = mu + v * lambda;
    mpz_class l = k * t - h;
    mpz_div(l.get_mpz_t(), l.get_mpz_t(), s.get_mpz_t());
    m = t * u * k - h * u - x.c * s;
    mpz_div(m.get_mpz_t(), m.get_mpz_t(), mpz_class (s * t).get_mpz_t());

    // A = st
    // B = ju - kt + ls
    // C = kl - jm
    a = s * t;
    b = j * u - k * t - l * s;
    mpz_class c = k * l - j * m;
    return elem({a, b, c});
}

classgroup::classElem classgroup::id_(mpz_class d) {
    mpz_class a = 1;
    mpz_class b = 1;

    // c = (b * b - d) / 4a
    mpz_class c = (1 - d) / 4;
    return { a, b, c };
}

classgroup::classElem classgroup::inv_(classgroup::classElem x) {
    return {x.a, -x.b, x.c};
}

classgroup::classElem classgroup::exp_(classgroup::classElem a, mpz_class n) {
    classgroup::classElem val;
    if(n < 0){
        val = classgroup::id();
        a = classgroup::inv(a);
        n = -n;
    } else {
        val = classgroup::id();
        a = a;
        n = n;
    }

    do{
        if(n == 0){
            return val;
        }
        if(mpz_odd_p(n.get_mpz_t())){
            val = classgroup::op(val, a);
        }
        a = classgroup::square(a);
        n >>= 1;
    } while (1);
}

classgroup::classElem classgroup::unknown_order_elem_(mpz_class d) {
    mpz_class a = 2;
    mpz_class b = 1;
    mpz_class c = (1 - d) / 8;
    return { a, b, c };
}

bool classgroup::eq(classgroup::classElem self, classgroup::classElem other) {
    return self.a == other.a && self.b == other.b && self.c == other.c;
}

classgroup::classElem classgroup::elem(classgroup::classElem abc) {
    classElem elem;
    elem = classgroup::reduce(abc.a, abc.b, abc.c);
    assert(validate(abc.a, abc.b, abc.c));
    return elem;
}

classgroup::classElem classgroup::construct_raw_elem_from_strings(char *a, char *b, char *c) {
    return {mpz_class (a), mpz_class (b), mpz_class (c)};
}

void classgroup::test_bad_elem() {
    printf("1\n");
    classgroup::elem({1, 2, 3});
}

void classgroup::test_elem_from() {
    mpz_class a1 = mpz_class ("16");
    mpz_class b1 = mpz_class ("105");
    mpz_class c1 = mpz_class (
            "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
      672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
      531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
      198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
      494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
      038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
      9057462766047140854869124473221137588347335081555186814207"
    );

    mpz_class a2 = mpz_class ("16");
    mpz_class b2 = mpz_class ("9");
    mpz_class c2 = mpz_class (
            "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
      672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
      531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
      198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
      494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
      038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
      9057462766047140854869124473221137588347335081555186814036"
    );

    classElem reduced_elem = classgroup::elem({a1, b1, c1});
    classElem also_reduced_elem = classgroup::elem({a2, b2, c2});
    assert(classgroup::eq(reduced_elem, also_reduced_elem));
}

void classgroup::test_equality() {
    classElem not_reduced = construct_raw_elem_from_strings(
            "16",
            "105",
            "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
      672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
      531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
      198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
      494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
      038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
      9057462766047140854869124473221137588347335081555186814207"
    );

    classElem reduced_ground_truth = construct_raw_elem_from_strings(
            "16",
            "9",
            "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
      672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
      531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
      198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
      494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
      038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
      9057462766047140854869124473221137588347335081555186814036"
    );

    classElem diff_elem = construct_raw_elem_from_strings(
            "4",
            "1",
            "19135043146754702466933535947700509509683047735103167439198117911126500023332446530136407244\
      268818886844950990589400824048541137030123517295048625578863052039394052606960429510076477727\
      812619793559333896857655440664448190570209733309248852860771133554929587999582285331513410741\
      679548532925359795754799072731031327175516868367484717873943724678975890638662600637655379895\
      797691446827331865910685793896910463236233398285859677535633644394859647446063344540995395360\
      815557919878168193309083573295900545539758915028677094752412489256178770608972743880695597825\
      16229851064188563419476497892884550353389340326220747256139"
    );

    assert(!eq(not_reduced, reduced_ground_truth));
    assert(!eq(not_reduced, diff_elem));
    assert(!eq(reduced_ground_truth, diff_elem));

    classElem reduced = classgroup::elem({not_reduced.a, not_reduced.b, not_reduced.c});
    assert(eq(reduced, reduced_ground_truth));
}

void classgroup::test_reduce_basic() {
    // Unreduced element.
    classElem to_reduce = construct_raw_elem_from_strings(
            "59162244921619725812008939143220718157267937427074598447911241410131470159247784852210767449\
      675610037288729551814191198624164179866076352187405442496568188988272422133088755036699145362\
      385840772236403043664778415471196678638241785773530531198720497580622741709880533724904220122\
      358854068046553219863419609777498761804625479650772123754523807001976654588225908928022367436\
      8",
            "18760351095004839755193532164856605650590306627169248964100884295652838905828158941233738613\
      175821849253748329102319504958410190952820220503570113920576542676928659211807590199941027958\
      195895385446372444261885022800653454209101497963588809819572703579484085278913354621371362285\
      341138299691587953249270188429393417132110841259813122945626515477865766896056280729710478647\
      13",
            "14872270891432803054791175727694631095755964943358394411314110783404577714102170379700365256\
      599679049493824862742803590079461712691146098397470840896560034332315858221821103076776907123\
      277315116632337385101204055232891361405428635972040596205450316747012080794838691280547894128\
      246741601088755087359234554141346980837292342320288111397175220296098629890108459305643419353\
      36"
    );

    classElem reduced_ground_truth = construct_raw_elem_from_strings(
            "26888935961824081232597112540509824504614070059776273347136888921115497522070287009841688662\
      983066376019079593372296556420848446780369918809384119124783870290778875424468497961559643807\
      918398860928578027038014112641529893817109240852544158309292025321122680747989987560029531021\
      808743313150630063377037854944",
            "14529985196481999393995154363327100184407232892559561136140792409262328867440167480822808496\
      853924547751298342980606034124112579835255733824790020119078588372593288210628255956605240171\
      744703418426092073347584357826862813733154338737148962212641444735717023402201569115323580814\
      54099903972209626147819759991",
            "28467266502267127591420289007165819749231433586093061478772560429058231137856046130384492811\
      816456933286039468940950129263300933723839212086399375780796041634531383342902918719073416087\
      614456845205980227091403964285870107268917183244016635907926846271829374679124848388403486656\
      1564478239095738726823372184204"
    );

    classElem res = classgroup::reduce(to_reduce.a, to_reduce.b, to_reduce.c);
    assert(eq(res, reduced_ground_truth));

    res = classgroup::reduce(reduced_ground_truth.a, reduced_ground_truth.b, reduced_ground_truth.c);
    assert(eq(res, reduced_ground_truth));
}

void classgroup::test_normalize_basic() {
    classElem unnormalized = construct_raw_elem_from_strings(
            "16",
            "105",
            "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
      672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
      531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
      198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
      494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
      038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
      9057462766047140854869124473221137588347335081555186814207"
    );

    classElem normalized_ground_truth = construct_raw_elem_from_strings(
            "16",
            "9",
            "4783760786688675616733383986925127377420761933775791859799529477781625005833111632534101811\
       06720472171123774764735020601213528425753087932376215639471576300984851315174010737751911943\
       19531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526\
       85419887133231339948938699768182757831793879217091871179468485931169743972659665650159413844\
       97394942286170683296647767144847422761580905834957146491938390841109871491186151583613524884\
       88402038894799695420483272708933239751363849397287571692736881031223140446926522431859701738\
       9945629057462766047140854869124473221137588347335081555186814036"
    );

    classElem res = classgroup::normalize(unnormalized.a, unnormalized.b, unnormalized.c);
    assert(eq(normalized_ground_truth, res));
}

void classgroup::test_discriminant_basic() {
    classElem g = classgroup::unknown_order_elem();
    assert(classgroup::discriminant(g.a, g.b, g.c) == DETAL);
}

void classgroup::test_discriminant_across_ops() {
    classElem id = classgroup::id();
    classElem g1 = classgroup::unknown_order_elem();
    classElem g2 = classgroup::op(g1, g1);
    classElem g3 = classgroup::op(id, g2);
    classElem g3_inv = classgroup::inv(g3);

    assert(classgroup::validate(id.a, id.b, id.c));
    assert(classgroup::validate(g1.a, g1.b, g1.c));
    assert(classgroup::validate(g2.a, g2.b, g2.c));
    assert(classgroup::validate(g3.a, g3.b, g3.c));
    assert(classgroup::validate(g3_inv.a, g3_inv.b, g3_inv.c));
}

void classgroup::test_op_single() {
    classElem a = construct_raw_elem_from_strings(
            "4",
            "1",
            "19135043146754702466933535947700509509683047735103167439198117911126500023332446530136407244\
      268818886844950990589400824048541137030123517295048625578863052039394052606960429510076477727\
      812619793559333896857655440664448190570209733309248852860771133554929587999582285331513410741\
      679548532925359795754799072731031327175516868367484717873943724678975890638662600637655379895\
      797691446827331865910685793896910463236233398285859677535633644394859647446063344540995395360\
      815557919878168193309083573295900545539758915028677094752412489256178770608972743880695597825\
      16229851064188563419476497892884550353389340326220747256139"
    );

    classElem b = construct_raw_elem_from_strings(
            "16",
            "41",
            "47837607866886756167333839869251273774207619337757918597995294777816250058331116325341018110\
      672047217112377476473502060121352842575308793237621563947157630098485131517401073775191194319\
      531549483898334742144138601661120476425524333273122132151927833887323969998955713328783526854\
      198871332313399489386997681827578317938792170918711794684859311697439726596656501594138449739\
      494228617068329664776714484742276158090583495714649193839084110987149118615158361352488488402\
      038894799695420483272708933239751363849397287571692736881031223140446926522431859701738994562\
      9057462766047140854869124473221137588347335081555186814061"
    );

    classElem ground_truth = construct_raw_elem_from_strings(
            "64",
            "9",
            "11959401966721689041833459967312818443551904834439479649498823694454062514582779081335254527\
      668011804278094369118375515030338210643827198309405390986789407524621282879350268443797798579\
      882887370974583685536034650415280119106381083318280533037981958471830992499738928332195881713\
      549717833078349872346749420456894579484698042729677948671214827924359931649164125398534612434\
      873557154267082416194178621185569039522645873928662298459771027746787279653789590338122122100\
      50972369992385512081817723330993784096234932189292318422025780578511173163060796492543474864\
      07264365691511785213717281118305284397086833770388796703509"
    );

    assert(eq(classgroup::op(a, b), ground_truth));
}

void classgroup::test_op_alternating() {
    classElem g_anchor = classgroup::unknown_order_elem();
    classElem g = classgroup::id();
    classElem g_star = classgroup::id();

    // g
    g = classgroup::op(g_anchor, g);

    // g^2, g^* = g^2
    g = classgroup::op(g_anchor, g);
    g_star = classgroup::op(g, g_star);

    // g^3
    g = classgroup::op(g_anchor, g);

    // g^4, g^* = g^2 * g^4 = g^6
    g = classgroup::op(g_anchor, g);
    g_star = classgroup::op(g, g_star);

    classElem ground_truth = construct_raw_elem_from_strings(
            "64",
            "9",
            "11959401966721689041833459967312818443551904834439479649498823694454062514582779081335254527\
      668011804278094369118375515030338210643827198309405390986789407524621282879350268443797798579\
      882887370974583685536034650415280119106381083318280533037981958471830992499738928332195881713\
      549717833078349872346749420456894579484698042729677948671214827924359931649164125398534612434\
      873557154267082416194178621185569039522645873928662298459771027746787279653789590338122122100\
      509723699923855120818177233309937840962349321892923184220257805785111731630607964925434748640\
      7264365691511785213717281118305284397086833770388796703509"
    );

    assert(eq(ground_truth, g_star));
}

void classgroup::test_op_complex() {
    // 1. Take g^100, g^200, ..., g^1000.
    // 2. Compute g^* = g^100 * ... * g^1000.
    // 3. For each of g^100, g^200, ..., g^1000 compute the inverse of that element and assert that
    //    g^* * current_inverse = product of g^100, g^200, ..., g^1000 without the inversed-out
    //    element.
    classElem g_anchor = classgroup::unknown_order_elem();
    classElem g = classgroup::id();

    vector<classElem> gs;
    vector<classElem> gs_invs;

    classElem g_star = classgroup::id();
    for(int i = 1; i <= 1000; i++){
        g = classgroup::op(g_anchor, g);
        assert(classgroup::validate(g.a, g.b, g.c));
        if (i % 100 == 0 ){
            gs.push_back(g);
            gs_invs.push_back(classgroup::inv(g));
            g_star = classgroup::op(g, g_star);
            assert(classgroup::validate(g_star.a, g_star.b, g_star.c));
        }
    }
    for(int i = 0; i < gs.size(); i++){
        assert(classgroup::validate(gs[i].a, gs[i].b, gs[i].c));
        assert(classgroup::validate(gs_invs[i].a, gs_invs[i].b, gs_invs[i].c));
        classElem curr_prod = classgroup::id();
        for (auto elem : gs) {
                if (!eq(elem, gs[i])) {
                    curr_prod = classgroup::op(curr_prod, elem);
                    assert(classgroup::validate(
                            curr_prod.a,
                            curr_prod.b,
                            curr_prod.c
                    ));
                }
        }
        assert(eq(classgroup::id(), classgroup::op(gs_invs[i], gs[i])));
        assert(eq(curr_prod, classgroup::op(gs_invs[i], g_star)));
    }

}

void classgroup::test_id_basic() {
    classElem g = classgroup::unknown_order_elem();
    classElem id = classgroup::id();
    assert(eq(g, classgroup::op(g, id)));
    assert(eq(g, classgroup::op(id, g)));
    assert(eq(id, classgroup::op(id, id)));
}

void classgroup::test_id_repeated() {
    classElem id = classgroup::id();
    classElem g_anchor = classgroup::unknown_order_elem();
    classElem g = classgroup::unknown_order_elem();
    for(int i = 0; i < 1000; i++) {
        id = classgroup::op(id, id);
        assert(eq(id, classgroup::id()));
        g = classgroup::op(g, classgroup::id());
        assert(eq(g, g_anchor));
    }
}

void classgroup::test_inv() {
    classElem id = classgroup::id();
    classElem g_anchor = classgroup::unknown_order_elem();
    classElem g = classgroup::unknown_order_elem();

    for(int i = 0; i < 1000; i++) {
        g = classgroup::op(g, g_anchor);
        classElem g_inv = classgroup::inv(g);
        assert(eq(id, classgroup::op(g_inv, g)));
        assert(eq(id, classgroup::op(g, g_inv)));
        assert(eq(g, classgroup::inv(g_inv)));
    }
}

void classgroup::test_exp_basic() {
    classElem g_anchor = classgroup::unknown_order_elem();
    classElem g = classgroup::id();

    for (int i = 1; i <= 1000; i++) {
        g = classgroup::op(g, g_anchor);
        assert(eq(g, classgroup::exp(g_anchor, i)));
    }
}

void classgroup::test_square_basic() {
    classElem g = classgroup::unknown_order_elem();
    classElem g4 = classgroup::id();

    // g^4
    for(int i = 0; i < 4; i++) {
        g4 = classgroup::op(g, g4);
    }

    // g^2
    classElem g2 = classgroup::op(g, g);

    // g^4
    g2 = classgroup::square(g2);

    assert(eq(g2, g4));
}







